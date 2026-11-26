# Node wrapper (design 014), sourced by every rule of an api >= 0.8.0
# Snakefile from the output root. `ob_run CMD...` runs a module and, when it
# fails, writes failure.json into the node directory. Bash and coreutils only:
# this runs inside the module's environment, which may have no Python.
#
# Set by the rule: OB_NODE (node dir), OB_UPSTREAM (producer dirs), OB_OUTPUTS
# (absolute output paths), OB_MANIFEST (a gather's lineage.json), LOG_FILE.
# Set by `ob run`: OB_TIMEOUT (seconds, 0 = none), OB_CONTINUE_ON_ERROR (1 = -k).
OB_OUT=$PWD

# Messages meaning the module never got going: missing script, package or
# library (Python, R, the dynamic linker).
OB_START_PATTERNS="ModuleNotFoundError|ImportError|can't open file|there is no package called|cannot open file|error while loading shared libraries"

_ob_json() {
    local s=${1//\\/\\\\}
    s=${s//\"/\\\"}
    printf '"%s"' "$(printf '%s' "$s" | tr -d '\000-\037' | cut -c1-300)"
}

# _ob_fail KIND EXIT_CODE REASON [UPSTREAM_DIRS]: record, then exit.
_ob_fail() {
    local kind=$1 rc=$2 reason=$3 d up="" t=${OB_TIMEOUT:-0} elapsed=$SECONDS
    for d in ${4:-}; do up+="${up:+, }$(_ob_json "$d")"; done
    [ -n "$rc" ] || elapsed=null
    printf '{"kind": "%s", "reason": %s, "exit_code": %s, "elapsed_s": %s, "timeout_s": %s, "upstream": [%s], "log": %s}\n' \
        "$kind" "$(_ob_json "$reason")" "${rc:-null}" "$elapsed" \
        "$([ "$t" -gt 0 ] && echo "$t" || echo null)" "$up" \
        "$(_ob_json "${LOG_FILE:+${LOG_FILE#"$OB_OUT"/}}")" > "$OB_OUT/$OB_NODE/failure.json"
    # Outputs exist for Snakemake, empty so nothing reads them as a result;
    # consumers go by failure.json, never by size.
    for d in ${OB_OUTPUTS:-}; do [ -d "$d" ] || : > "$d"; done
    echo "ob: $kind: $reason" >&2
    [ "${OB_CONTINUE_ON_ERROR:-0}" = 1 ] && exit 0
    exit $(( ${rc:-1} > 0 ? ${rc:-1} : 1 ))
}

# _ob_timeout SECONDS CMD...: GNU `timeout -k 30`, which macOS lacks. Exits
# 124 when the time ran out (KILL 30 s after TERM), else with CMD's status.
# A subshell, so `set -m` (CMD in its own process group, signalled as a
# whole like `timeout` does) stays local.
_ob_timeout() (
    t=$1
    shift
    set -m +e
    flag=$(mktemp)
    rm -f "$flag"
    "$@" &
    pid=$!
    # The watcher's sleeps are killed with it, so none outlives the module.
    (
        trap 'kill $s 2>/dev/null; exit' TERM
        sleep "$t" & s=$!; wait $s
        : > "$flag"
        kill -TERM -- -$pid
        sleep 30 & s=$!; wait $s
        kill -KILL -- -$pid
    ) > /dev/null 2>&1 &
    watcher=$!
    wait $pid 2> /dev/null  # the job notice ("Terminated"), not CMD's stderr
    rc=$?
    kill $watcher 2> /dev/null
    wait $watcher 2> /dev/null
    if [ -e "$flag" ]; then
        rm -f "$flag"
        exit 124
    fi
    exit $rc
)

ob_run() {
    local d rc up="" t=${OB_TIMEOUT:-0} err fifo teepid msg cmd=("$@")
    rm -f "$OB_OUT/$OB_NODE/failure.json"
    for d in ${OB_UPSTREAM:-}; do
        [ -e "$OB_OUT/$d/failure.json" ] && up+=" $d"
    done
    [ -n "$up" ] && _ob_fail upstream "" "input from a failed node" "$up"
    if [ -n "${OB_MANIFEST:-}" ] && grep -q '"members": \[\]' "$OB_MANIFEST"; then
        _ob_fail upstream "" "every member failed"
    fi

    [ "$t" -gt 0 ] && cmd=(_ob_timeout "$t" "$@")
    err=$(mktemp)
    SECONDS=0
    set +e
    # stderr to the log as usual, and to $err for the reason. Through a fifo,
    # not a pipe: a pipe waits for every process holding its write end, so a
    # descendant that outlives the module (a daemon, or one stuck in the
    # kernel where even KILL cannot reach it) would hold the node open past
    # its timeout. tee gets 5 s to drain once the module has returned.
    fifo="$err.fifo"
    mkfifo "$fifo"
    tee "$err" < "$fifo" >&2 &
    teepid=$!
    "${cmd[@]}" 2> "$fifo"
    rc=$?
    for _ in 1 2 3 4 5 6 7 8 9 10; do
        kill -0 "$teepid" 2> /dev/null || break
        sleep 0.5
    done
    kill "$teepid" 2> /dev/null
    wait "$teepid" 2> /dev/null
    rm -f "$fifo"
    set -e
    msg=$(grep -E -m1 "$OB_START_PATTERNS" "$err" || true)
    [ -n "$msg" ] || msg=$(grep -v -e '^[[:space:]]*$' -e '^Execution halted' "$err" | tail -n1 || true)
    rm -f "$err"

    if [ "$t" -gt 0 ] && [ "$rc" = 124 ]; then
        _ob_fail timed_out "$rc" "timed out after $t s"
    elif [ "$rc" = 126 ]; then
        _ob_fail failed_to_start "$rc" "entrypoint not executable"
    elif [ "$rc" = 127 ]; then
        _ob_fail failed_to_start "$rc" "entrypoint or interpreter not found"
    elif [ "$rc" != 0 ] && grep -qE "$OB_START_PATTERNS" <<< "$msg"; then
        _ob_fail failed_to_start "$rc" "$msg"
    elif [ "$rc" -gt 128 ]; then
        _ob_fail crashed "$rc" "killed by signal $((rc - 128))$([ "$rc" = 137 ] && echo ' (out of memory?)')"
    elif [ "$rc" != 0 ]; then
        _ob_fail crashed "$rc" "exit code $rc${msg:+: $msg}"
    fi
    for d in ${OB_OUTPUTS:-}; do
        [ -e "$d" ] || _ob_fail missing_output 0 "declared output $(basename "$d") not written"
    done
}
