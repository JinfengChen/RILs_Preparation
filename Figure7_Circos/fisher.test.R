print("Somatic VS. Strains")
x <- matrix(c(292, 99, 56, 10),  nrow = 2)
fisher.test(x)
chisq.test(x)

print("RIL VS. Strains")
x <- matrix(c(723, 99, 165, 10),  nrow = 2)
fisher.test(x)
chisq.test(x)

print("Somatic VS. RIL")
x <- matrix(c(292, 723, 56, 165),  nrow = 2)
fisher.test(x)
chisq.test(x)

print("Unique VS. Simulation")
x <- matrix(c(881, 130, 186, 141),  nrow = 2)
fisher.test(x)
chisq.test(x)

