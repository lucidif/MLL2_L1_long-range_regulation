process JUICER_PRE {
    tag "${meta.id}"
    label 'process_medium'
    
    container "'lucidif/microc:0.0.1'"

    input:
    tuple val(meta), path(pairs)
    path hic_tools_jar
    path chromsize
    val res

    output:
    tuple val(meta), path("*.hic")               , emit: hic
    path "versions.yml"                          , emit: versions

    script:
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def avail_mem = 4
    if (!task.memory) {
        log.info 'Available memory not known - defaulting to 4GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    java -Xms512m -Xmx${avail_mem}g \\
        -jar ${hic_tools_jar} pre \\
        -r ${res} \\
        $args \\
        --threads $task.cpus \\
        ${pairs} ${prefix}.${res}.hic ${chromsize}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(echo \$(java -jar ${hic_tools_jar} --version 2>&1) | sed 's/^.*Version //; s/Usage.*\$//')
    END_VERSIONS
    """
}
