set weights {}
set biases {}
set reactant_center [list -75.0 75.0]
set product_center [list 80.0 -80.0]
set radius 10.0

proc NC_init {param_file} {
    global weights biases

    source $param_file
    foreach var {fc_0_weight fc_2_weight fc_4_weight fc_6_weight fc_0_bias fc_2_bias fc_4_bias fc_6_bias} {
        if {![info exists $var]} {
            error "Missing parameter: $var"
        }
    }

    set weights [list $fc_0_weight $fc_2_weight $fc_4_weight $fc_6_weight]
    set biases [list $fc_0_bias $fc_2_bias $fc_4_bias $fc_6_bias]
}

proc _smooth_transition {x center radius} {
    set d2 0.0
    foreach xi $x ci $center {
        set xi_clean [regsub -all {[^0-9.eE+-]} $xi ""]
        set ci_clean [regsub -all {[^0-9.eE+-]} $ci ""]

        set diff [expr {double($xi_clean) - double($ci_clean)}]
        set d2 [expr {$d2 + $diff**2}]
    }
    set r_plus_sq [expr {(double($radius) + 2) ** 2}]
    set arg [expr {1000.0 * ($d2 - $r_plus_sq)}]
    return [expr {0.5 - 0.5*tanh($arg)}]
}

proc _d_smooth_transition {x center radius} {
    set dx {}
    foreach xi $x ci $center {
        set xi_clean [regsub -all {[^0-9.eE+-]} $xi ""]
        set ci_clean [regsub -all {[^0-9.eE+-]} $ci ""]

        lappend dx [expr {double($xi_clean) - double($ci_clean)}]
    }
    set d2 [expr {[lindex $dx 0] ** 2 + [lindex $dx 1] ** 2}]
    set r_plus_sq [expr {(double($radius) + 2) ** 2}]
    set arg [expr {1000.0 * ($d2 - $r_plus_sq)}]
    set sech2 [expr {1.0 / (cosh($arg)**2)}]
    set result {}
    foreach d $dx {
        lappend result [expr {-1000.0 * $d * $sech2}]
    }
    return $result
}

proc NC_sigmoid {x} {
    return [expr {1.0 / (1.0 + exp(-$x))}]
}

proc NC_vector_add {v1 v2} {
    set result {}
    foreach a $v1 b $v2 {
        lappend result [expr {$a + $b}]
    }
    return $result
}

proc NC_vector_matrix_mult {vec matrix} {
    set result {}
    foreach row $matrix {
        if {[llength $vec] != [llength $row]} {
            error "Dimension mismatch in vector_matrix_mult"
        }

        set sum 0.0
        foreach v $vec r $row {
            set sum [expr {$sum + $v * $r}]
        }
        lappend result $sum
    }
    return $result
}

proc NC_transpose {matrix} {
    set transposed {}
    set cols [llength [lindex $matrix 0]]
    for {set i 0} {$i < $cols} {incr i} {
        set new_row {}
        foreach row $matrix {
            lappend new_row [lindex $row $i]
        }
        lappend transposed $new_row
    }
    return $transposed
}

proc NC_forward {x} {
    global weights biases reactant_center product_center radius
    set rhoA [_smooth_transition $x $reactant_center $radius]
    set rhoB [_smooth_transition $x $product_center $radius]
    #set h $x
    set x_nor {}
    foreach val $x {
        lappend x_nor [expr {double($val)/180.0}]
    }
    set h $x_nor
    for {set i 0} {$i < 3} {incr i} {
        set h [NC_vector_matrix_mult $h [lindex $weights $i]]
        set h [NC_vector_add $h [lindex $biases $i]]
        #set h [lmap val $h {expr {tanh($val)}}]
        set h_temp {}
        foreach val $h {
            set result [expr {tanh($val)}]
            lappend h_temp $result
        }
        set h $h_temp
    }
    set h [NC_vector_matrix_mult $h [lindex $weights 3]]
    set h [NC_vector_add $h [lindex $biases 3]]
    set u [NC_sigmoid [lindex $h 0]]
    set q [expr {(1.0 - $rhoA - $rhoB) * $u + $rhoB}]
    return $q
}


#NC_init "network_params.tcl"
#set position [list -75 0]
#set q_value [NC_forward $position]
#set gradient [NC_gradient $position]
#puts "Committor value: [format "%.6e" $q_value]"
#puts "Gradient: $gradient"
