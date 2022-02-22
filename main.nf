params.aaseqs = "testdata/seqs.pep"
params.outdir = "results"

log.info """\
            _______  _______                    _
  |\     /|(  ___  )(  ____ )                  ( )
  | )   ( || (   ) || (    )|                  | |
  | (___) || (___) || (____)|                  | |
  |  ___  ||  ___  ||  _____) _______          | |
  | (   ) || (   ) || (      (  ____ )|\     /|(_)
  | )   ( || )   ( || )  _   | (    )|( \   / ) _
  |/     \||/     \||/  (_)  | (____)| \ (_) / (_)
                             |  _____)  \   /
                             | (         ) (
                             | )         | |
                             |/          \_/
Running happy with the options:
aaseqs = ${params.aaseqs}
outdir = ${params.outdir}
"""