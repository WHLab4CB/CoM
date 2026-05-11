#!/usr/bin/env tclsh

### read the coordinates of anchors
set ank_id_list [list]
set ank_phi_list [list]
set ank_psi_list [list]
set f_ank [open "./ank.txt" r]
while {[gets $f_ank line] != -1} {
    set line_data [split $line " "]  
    set ank_id [lindex $line_data 0]    
    set ank_phi [lindex $line_data 1]
    set ank_psi [lindex $line_data 2]
    lappend ank_id_list $ank_id
    lappend ank_phi_list $ank_phi
    lappend ank_psi_list $ank_psi
}
close $f_ank
set n_ank [llength $ank_id_list]

### params
set reactant_phi -75 
set reactant_psi 75
set product_phi 80
set product_psi -80
set a_radius 10
set ini_ank_id 32
set k_con 1


### calculate bias force
proc cal_bias_force {colvar_phi colvar_psi ini_ank_id ank_phi_list ank_psi_list n_ank k_con} {
    set this_phi [lindex $ank_phi_list $ini_ank_id]
    set this_psi [lindex $ank_psi_list $ini_ank_id]
    
    set dist_this_1 [expr {$colvar_phi-$this_phi}]
    set dist_this_2 [expr {$colvar_psi-$this_psi}]
    set dist0 [expr {sqrt($dist_this_1*$dist_this_1 + $dist_this_2*$dist_this_2)}]

    set bf_phi 0.0
    set bf_psi 0.0
    for {set i 0} {$i < $n_ank} {incr i} {
        if {$i == $ini_ank_id} {
            continue
        } else {
            set dist_this_3 [expr {$colvar_phi-[lindex $ank_phi_list $i]}]
            set dist_this_4 [expr {$colvar_psi-[lindex $ank_psi_list $i]}]
            set dist1 [expr {sqrt($dist_this_3*$dist_this_3 + $dist_this_4*$dist_this_4)}]
            set dist_diff [expr {$dist0 - $dist1}]
            if {$dist_diff > 0} {
                set bf_phi [expr {$bf_phi-2.0*$k_con*$dist_diff*($dist_this_1/$dist0 - $dist_this_3/$dist1)}]
                set bf_psi [expr {$bf_psi-2.0*$k_con*$dist_diff*($dist_this_2/$dist0 - $dist_this_4/$dist1)}]
            }        
        }
    }
    return [list $bf_phi $bf_psi]
}


proc calc_colvar_forces {t} {
    global reactant_phi reactant_psi product_phi product_psi a_radius ini_ank_id k_con ank_id_list ank_phi_list ank_psi_list n_ank

    if {$t > 0} {
        set cov_phi [cv colvar phi value]
        set cov_psi [cv colvar psi value]
        
        set bf_all [cal_bias_force $cov_phi $cov_psi $ini_ank_id $ank_phi_list $ank_psi_list $n_ank $k_con]
        set bf_phi [lindex $bf_all 0]
        set bf_psi [lindex $bf_all 1]

        cv colvar phi addforce $bf_phi
        cv colvar psi addforce $bf_psi
    }
}



 
