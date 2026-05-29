source "./committor.tcl"
NC_init "./network_params.tcl"
set qlist [list 0 3e-5 3e-4 3e-3 3e-2 0.1 0.3 0.5 0.6 0.7 0.8 0.9 1.0]
set qthis 4
set qlen [llength $qlist]
set qtarg [lindex $qlist $qthis]
if {$qthis == 0} {
    set qleft -1
    set qright [expr {$qthis + 1}]
} elseif {$qthis == [expr {$qlen - 1}]} {
    set qleft [expr {$qthis - 1}]
    set qright -2
} else {
    set qleft [expr {$qthis - 1}]
    set qright [expr {$qthis + 1}]
}

proc calcforces {} {
    global qlist qthis qlen qtarg qleft qright reactant_center product_center radius

    set t [getstep]
    if {$t>0} { 
        set cov_h1 [cv colvar hbond1 value]
        set cov_h2 [cv colvar hbond2 value]
        set cov_h3 [cv colvar hbond3 value]
        set cov_h4 [cv colvar hbond4 value]
        set cov_sum_h [cv colvar dist_hbond value]
        set q_input [list $cov_h1 $cov_h2 $cov_h3 $cov_h4]
        set q_value [NC_forward $q_input]

        if {$qthis == 0} {
            if {$q_value > [lindex $qlist $qright]} {
                set t_final [expr $t/1000.0]
                print "final_result $qthis $qright $t_final $q_value"
                exit
            }
        } elseif {$qthis == [expr {$qlen - 1}]} {
            if {$q_value < [lindex $qlist $qleft]} {
                set t_final [expr $t/1000.0]
                print "final_result $qthis $qleft $t_final $q_value"
                exit
            }
        } else {
            if {$qleft == 0} {
                if {$cov_sum_h < $reactant_center} {
                    set t_final [expr $t/1000.0]
                    print "final_result $qthis $qleft $t_final $q_value"
                    exit
                }
            } else {
                if {$q_value < [lindex $qlist $qleft]} {
                    set t_final [expr $t/1000.0]
                    print "final_result $qthis $qleft $t_final $q_value"
                    exit
                }
            }
            if {$qright == [expr {$qlen - 1}]} {
                if {$cov_sum_h > $product_center} {
                    set t_final [expr $t/1000.0]
                    print "final_result $qthis $qright $t_final $q_value"
                    exit
                }
            } else {
                if {$q_value > [lindex $qlist $qright]} {
                    set t_final [expr $t/1000.0]
                    print "final_result $qthis $qright $t_final $q_value"
                    exit
                }
            }
        }
    }
}
