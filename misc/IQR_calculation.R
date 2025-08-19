library(gridExtra)
library(ggplot2)

hg19_hg18=read.delim("/sc/arion/projects/mscic1/results/anina/fun_project_4.23/Rosmap_combos/hg19_hg18/DE/DE_assembly_corrected_for_id_assembly_test.txt")
v30_hg18=read.delim("/sc/arion/projects/mscic1/results/anina/fun_project_4.23/Rosmap_combos/v30_hg18/DE/DE_assembly_corrected_for_id_assembly_test.txt")
v30_hg19=read.delim("/sc/arion/projects/mscic1/results/anina/fun_project_4.23/Rosmap_combos/v30_hg19/DE/DE_assembly_corrected_for_id_assembly_test.txt")
v30_v43=read.delim("/sc/arion/projects/mscic1/results/anina/fun_project_4.23/Rosmap_combos/v30_v43/DE/DE_assembly_corrected_for_id_assembly_test.txt")
v43_hg18=read.delim("/sc/arion/projects/mscic1/results/anina/fun_project_4.23/Rosmap_combos/v43_hg18/DE/DE_assembly_corrected_for_id_assembly_test.txt")
v43_hg19=read.delim("/sc/arion/projects/mscic1/results/anina/fun_project_4.23/Rosmap_combos/v43_hg19/DE/DE_assembly_corrected_for_id_assembly_test.txt")
T2T_hg18=read.delim("/sc/arion/projects/mscic1/results/anina/fun_project_4.23/Rosmap_combos/T2T_hg18/DE/DE_assembly_corrected_for_id_assembly_test.txt")
T2T_hg19=read.delim("/sc/arion/projects/mscic1/results/anina/fun_project_4.23/Rosmap_combos/T2T_hg19/DE/DE_assembly_corrected_for_id_assembly_test.txt")
T2T_v30=read.delim("/sc/arion/projects/mscic1/results/anina/fun_project_4.23/Rosmap_combos/T2T_v30/DE/DE_assembly_corrected_for_id_assembly_test.txt")
T2T_v43=read.delim("/sc/arion/projects/mscic1/results/anina/fun_project_4.23/Rosmap_combos/T2T_v43/DE/DE_assembly_corrected_for_id_assembly_test.txt")

IQR(v30_v43$logFC)
# [1] 0.0008112102
IQR(T2T_v43$logFC)
# [1] 0.001830775
IQR(hg19_hg18$logFC)
# [1] 0.8765176
IQR(v30_hg18$logFC)
# [1] 2.09593
IQR(v30_hg19$logFC)
# [1] 0.01563655
IQR(v30_v43$logFC)
# [1] 0.0008112102
IQR(v43_hg18$logFC)
# [1] 5.185085
IQR(v43_hg19$logFC)
# [1] 0.1597644
IQR(T2T_hg18$logFC)
# [1] 2.968242
IQR(T2T_hg19$logFC)
# [1] 0.09545829
IQR(T2T_v30$logFC)
# [1] 0.002290513
IQR(T2T_v43$logFC)
# [1] 0.001830775

mean(c(IQR(v30_v43$logFC),IQR(T2T_v43$logFC),IQR(hg19_hg18$logFC),IQR(v30_hg18$logFC),IQR(v30_hg19$logFC),IQR(v30_v43$logFC),IQR(v43_hg18$logFC),IQR(v43_hg19$logFC),IQR(T2T_hg18$logFC),IQR(T2T_hg19$logFC),IQR(T2T_v30$logFC),IQR(T2T_v43$logFC) ))
#0.9503507
