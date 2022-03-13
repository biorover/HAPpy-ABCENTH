params.aaseqs = "$baseDir/testdata/seqs.pep"
params.outdir = "results"
params.cluster_threshold = 0.5
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

aaseqs = Channel.fromPath( params.aaseqs, checkIfExists: true )

process mafftree {
    
    input:
    path aaseqs 
     
    output:
    path 'mafft_tree.tre' 

    script:       
    """
    mafft --thread $task.cpus --treeout ${aaseqs}
    sed -e 's/[0-9]*_//' ${aaseqs}.tree > mafft_tree.tre
    """
}

process build_clusters {

    input:
    path tree
    path aaseqs

    output:
    path 'clusters'

    script:
    """
    mkdir clusters
    MAGOT build-clusters ${tree} ${aaseqs} clusters ${params.cluster_threshold}
    """
}


workflow {
    mafft_tree = mafftree(aaseqs)
    cluster_dir = build_clusters(mafft_tree,aaseqs)

}