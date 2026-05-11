#!/usr/bin/env tclsh


### params
set reactant_dist 15.0
set product_dist 40.0
set n_steps 3000

proc calcforces {} {
    global reactant_dist product_dist n_steps

    set t [getstep]
    if {$t>0} {
        set cov_dist [cv colvar dist_hbond value]

        ### exit and print
        if {$cov_dist < $reactant_dist || $cov_dist > $product_dist} {
            set cov_hbond1 [cv colvar hbond1 value]
            set cov_hbond2 [cv colvar hbond2 value]
            set cov_hbond3 [cv colvar hbond3 value]
            set cov_hbond4 [cv colvar hbond4 value]
            
            print "final_result $cov_dist $t $cov_hbond1 $cov_hbond2 $cov_hbond3 $cov_hbond4"
            exit
        }

        if {$t == $n_steps || $t == 1} {
            set cov_hbond1 [cv colvar hbond1 value]
            set cov_hbond2 [cv colvar hbond2 value]
            set cov_hbond3 [cv colvar hbond3 value]
            set cov_hbond4 [cv colvar hbond4 value]

            print "final_result $cov_dist $t $cov_hbond1 $cov_hbond2 $cov_hbond3 $cov_hbond4"
        }
    }    
}



 
