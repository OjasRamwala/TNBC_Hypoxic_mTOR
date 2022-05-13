files = list.files(pattern = ".bam.count")
count_matrix = data.frame(matrix(nrow = 64256))

for (i in files) {
  file = read.delim(i, row.names=1)
  count_matrix = cbind(count_matrix, file)
}
count_matrix = count_matrix[-c(1)]
count_matrix
