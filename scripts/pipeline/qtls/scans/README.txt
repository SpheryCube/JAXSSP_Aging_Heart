scan1(addcovar = sex, intcovar = sex)
scan1(addcovar = age, intcovar = age)

scan1(addcovar = sex)
scan1(addcovar = c(sex, age))

# scan1(addcovar = c(sex, age), intcovar = c(sex, age))



The qtl_submission_infrastructure folder contains files that help
submit qtl jobs to the HPC in parallel. Move those scripts out of that folder into this folder level before running the scans and submitting the jobs to the HPC.