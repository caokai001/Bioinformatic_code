'''
-*- coding: utf-8 -*-
@Author  : kcao
@Time    : 2019/5/29 10:19
@Software: PyCharm
@File    : script.py
'''
snp=r"C:\Users\16926\Desktop\项目\前列腺癌\snp_chip_ctcf\snp_indel\LNCaP-total-snp.bed"
fimo_gff=r"C:\Users\16926\Desktop\项目\前列腺癌\snp_chip_ctcf\Motif_fimo\LNCaP_CTCF_narrowfimo.gff\fimo.sort.gff"
peak=r"C:\Users\16926\Desktop\项目\前列腺癌\snp_chip_ctcf\CTCF_peak\LNCaP_CTCF.sort_narrow_peaks.narrowPeak"

snp_file=open(snp)
gff_file=open(fimo_gff)

new_file = open('LNCap.all.txt', 'w')


def fetch_line(file):
    line = file.readline()
    if not line:
        return None
    while line[0] == '#':
        line = file.readline()
    return line.rstrip().split('\t')


lines = None
with open(peak) as peaks:
    motif = None
    snp = None
    for peak in peaks:
        peak = peak.rstrip().split('\t')
        chrom = peak[0]
        st = int(peak[1])
        ed = int(peak[2])

        motifs = []  ###该peak里面motif 个数
        while True:
            if motif is None:
                motif = fetch_line(gff_file)
                if not motif:
                    break
            if motif[0] != chrom or int(motif[3]) >= ed:
                break
            else:
                motifs.append(';'.join([str(i) for i in motif]))
                print(motif[0], motif[3], motif[4])
                motif = None

        snps = []
        while True:
            if snp is None:
                snp = fetch_line(snp_file)
                if not snp:
                    break
            snp_pos = int(snp[1])
            if snp[0] != chrom:
                break
            if st <= snp_pos <= ed:
                snps.append(';'.join([str(i) for i in snp[:5]]))
                print(snp[0],snp[1],snp[2])
                snp = None
            elif snp_pos < st:
                snp = None
            else:
                break

        new_file.write("{}\t{}\t{}\n".format("\t".join(peak), len(motifs), len(snps)))
#
