params.aaseqs = "$baseDir/testdata/seqs.pep"
params.outdir = "results"
nextflow.enable.dsl=2

log.info """\
                _______  _______                    _
    |\\     /|(  ___  )(  ____ )                  ( )
    | )   ( || (   ) || (    )|                  | |
    | (___) || (___) || (____)|                  | |
    |  ___  ||  ___  ||  _____) _______          | |
    | (   ) || (   ) || (      (  ____ )|\\     /|(_)
    | )   ( || )   ( || )  _   | (    )|( \\   / ) _
    |/     \\||/     \\||/  (_)  | (____)| \\ (_) / (_)
                                |  _____)  \\   /
                                | (         ) (
                                | )         | |
                                |/          \\_/
    You are running happy with the following parameters:
    aaseqs = ${params.aaseqs}
    outdir = ${params.outdir}
    """.stripIndent()

process mafftree {
    
    input:
    path aaseqs from params.aaseqs
     
    output:
    path 'tree' into mafft_tree

    script:       
    """
    mafft --thread $task.cpus --treeout ${aaseqs}
    cp ${aaseqs}.tree tree
    """
}