source "./committor.tcl"
NC_init "./network_params.tcl"
set qlist [list 0 3e-4 1e-3 3e-3 1e-2 3e-2 0.1 0.5 0.9 0.98 0.995 1]
set qthis 0
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
        set cov_phi [cv colvar ang_phi value]
        set cov_psi [cv colvar ang_psi value]
        set q_input [list $cov_phi $cov_psi]
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
                set diff1 [expr {$cov_phi - [lindex $reactant_center 0]}]
                set diff2 [expr {$cov_psi - [lindex $reactant_center 1]}]
                set dist [expr {sqrt($diff1**2 + $diff2**2)}]
                if {$dist < $radius} {
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
                set diff1 [expr {$cov_phi - [lindex $product_center 0]}]
                set diff2 [expr {$cov_psi - [lindex $product_center 1]}]
                set dist [expr {sqrt($diff1**2 + $diff2**2)}]
                if {$dist < $radius} {
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
