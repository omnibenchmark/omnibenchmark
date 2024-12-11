-- R_4.4.1_Clustering

help([[
    R 4.4.1 environment environment with clustering libraries.
    Includes argparse, readr, dplyr, mclust, caret, and Python 3.12.6
]])

whatis("Name: R_4.4.1_Clustering")
whatis("Version: 1.1.0")
whatis("Category: OmniBenchmark")
whatis("Description: R environment for data clustering and analysis")

-- R Configuration
setenv("R_VERSION", "4.4.1")

-- Define the R environment path
local r_version = "4.4.1"
local r_home = capture("Rscript -e 'cat(Sys.getenv(\"R_HOME\"))'")
local r_bin = r_home .. "/bin"   -- R binaries are typically under R_HOME/bin
local r_lib = r_home .. "/library"  -- R library path is under R_HOME/library

-- Set up the environment to use R 4.4.1
prepend_path("PATH", r_bin)
prepend_path("LD_LIBRARY_PATH", r_lib)
prepend_path("R_LIBS", r_lib)

-- Load necessary R libraries (assuming they are installed in the specified path)
local required_libraries = {
    "argparse",
    "readr",
    "dplyr",
    "mclust",
    "caret"
}

-- Make sure the required R packages are accessible by adding their library path
for _, library in ipairs(required_libraries) do
    local library_path = r_lib .. "/" .. library
    prepend_path("R_LIBS", library_path)
end

-- Python Configuration
setenv("PYTHON_VERSION", "3.12.6")

-- Define the Python environment path
local python_version = "3.12"

local python_bin = capture("python" .. python_version .. " -c 'import sys; print(sys.executable)'")
local python_home = capture("python" .. python_version .. " -c 'import sys; print(sys.prefix)'")
local python_lib = python_home .. "/lib/python" .. python_version
local python_site_packages = python_lib .. "/site-packages"

-- Set up the environment to use Python 3.12.6
prepend_path("PATH", python_bin)
prepend_path("LD_LIBRARY_PATH", python_lib)
prepend_path("PYTHONPATH", python_site_packages)

-- Print information about the environment
LmodMessage("R 4.4.1 and clustering libraries (argparse, readr, dplyr, mclust, caret) along with Python 3.12.6 have been successfully loaded.")

---
-- Function to retrieve console output
--
function capture(cmd, raw)
    local handle = assert(io.popen(cmd, 'r'))
    local output = assert(handle:read('*a'))

    handle:close()

    if raw then
        return output
    end

    output = string.gsub(
        string.gsub(
            string.gsub(output, '^%s+', ''),
            '%s+$',
            ''
        ),
        '[\n\r]+',
        ' '
    )

   return output
end