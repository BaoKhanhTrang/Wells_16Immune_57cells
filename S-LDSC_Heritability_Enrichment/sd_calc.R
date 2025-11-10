input <- readLines(file("stdin"))
sd_value <- sd(input)
cat(sd_value,"\n")