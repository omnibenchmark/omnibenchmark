
def create_all_rule(paths):
    rule all:
        name: 'all'
        input:
            paths,
            # "out/data/D2/default/D2.data.ext",
            # "out/data/D2/default/process/P1/default/D2.txt.gz",
            # "out/data/D1/default/process/P2/default/methods/M2/default/D1.model.out.gz"
            # "out/data/D1/default/process/P2/default/methods/M2/default/m1/default/D1.results.txt"