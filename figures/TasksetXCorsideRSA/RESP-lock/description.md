Note - 5.22.26

From now on, want to be more conscientious abt naming and specifying what is what, esp since lots of diff versions of decoding etc. 

Files/plots should take the format:

modelV_(merged or not?)_bal_balanceV_evLock_extra.extension

where,
modelV is decoding variable
merged means data from all 3 sessions were concatenated for decoding
balanceV is what cross-validation folds were balanced on. so far, I have run decoding with:
 - NUMTASK
 - modelV (balanced by whatever modelV was
 - 'balance' - a newer variable that is each combination of TASKSET_CORSIDE x NUMTASK. 