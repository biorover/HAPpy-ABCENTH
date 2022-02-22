value_ch = Channel.value(3.14)
queue_ch = Channel.of( 1, 3, 5, 7 )
process basicExample {
    input:
    val x from value_ch
    val y from queue_ch
    shell:
    "echo $x $y"
}
