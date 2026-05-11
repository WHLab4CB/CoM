#!/usr/bin/env tclsh


### params
set reactant_phi -75 
set reactant_psi 75
set product_phi 80
set product_psi -80
set a_radius 10
set n_step 300

### calculate distance
proc cal_dist {colvar_phi colvar_psi center_phi center_psi} {
    set dist_phi [expr {$colvar_phi-$center_phi}]
    set dist_psi [expr {$colvar_psi-$center_psi}]
    set dist_0 [expr {$dist_phi*$dist_phi+$dist_psi*$dist_psi}]
    return $dist_0
}


proc calcforces {} {
    global reactant_phi reactant_psi product_phi product_psi a_radius n_step

    set t [getstep]
    if {$t>0} {
        set cov_phi [cv colvar phi value]
        set cov_psi [cv colvar psi value]
        
        set r_dist [cal_dist $cov_phi $cov_psi $reactant_phi $reactant_psi]
        set p_dist [cal_dist $cov_phi $cov_psi $product_phi $product_psi]

        if {$r_dist < $a_radius || $p_dist < $a_radius} {
            print "final_result $cov_phi $cov_psi $t $r_dist $p_dist"
            exit
        }

        if {$t == $n_step || $t == 1} {
            print "final_result $cov_phi $cov_psi $t $r_dist $p_dist"
        }


    }    
}



 
