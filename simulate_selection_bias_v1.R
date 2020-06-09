#*****************************************#
#       Test the selection bias           #
#*****************************************#

args = commandArgs(TRUE)
i = as.numeric(args[1])         # the i-th replicate
bxy = as.numeric(args[2])       # bxy of ACE on SCZ
f = as.numeric(args[3])         # allele frequency
gamma_x = as.numeric(args[4])   # effect of ACE on probability of selection bias
gamma_u = as.numeric(args[5])   # effect of confounding factor on probability of selection bias
rstfile = args[6]               # file to save results

#print(i)
#print(bxy)
#print(gamma_x)
#print(gamma_u)

# i = 1; # The ith simulation
# bxy = -0.55; 
# gamma_x = -0.2;
# gamma_u = -0.2;
# rstfile = "simu_slct_bias_bxy_-0.55_gamma_x_-0.2_gamma_u_-0.2_rep_1.txt"

# SMR analysis
SMR = function(bzx, bzx_se, bzx_pval, bzy, bzy_se) {
    bxy = bxy_se = bxy_pval = NA
    if(bzx_pval < 5e-8) {
        bxy = bzy/bzx
        bxy_se = sqrt( (bzy_se^2*bzx^2 + bzx_se^2*bzy^2) / bzx^4 )
        bxy_pval = pchisq((bxy/bxy_se)^2, 1, lower.tail=F)
    }
    return(list(bxy=bxy,bxy_se=bxy_se,bxy_pval=bxy_pval))
}

# set parameters
set.seed(i + f*100);
n1 = 32000; n2 = 500000  # sample sizes

hzx = 0.0025; bzx = sqrt(hzx/(2*f*(1-f)))  # h2 of x explained by z and bzx
bux = 1; hux = 0.3;     # h2 of x explained by u and bux, u - confounding factor
hxy = bxy^2             # h2 of y explained by x
buy = -1; huy = 0.3;    # h2 of y explained by u

gamma_0 = -1.6;         # logit(pai) = gamma_0 + gamma_x + gamma_u + e

p = 0.01; # population prevalence of SCZ

# dataset 1, eQTL, z -> ACE
z1 = rbinom(n1, 2, f);
u = scale(rnorm(n1, 0, 1))[,1]
x1 = z1*bzx + scale(bux*u)[,1]*sqrt(hux) + scale(rnorm(n1, 0, 1))[,1]*sqrt(1-hzx-hux); # x = z + u + e;
x1 = scale(x1)[,1]
# dataset 2, GWAS, z -> ACE -> SCZ
z2 = rbinom(n2, 2, f);
u = scale(rnorm(n2, 0, 1))[,1]
x2 = z2*bzx + scale(bux*u)[,1]*sqrt(hux) + scale(rnorm(n2, 0, 1))[,1]*sqrt(1-hzx-hux); # x = z + u + e;
x2 = scale(x2)[,1];
cov_xy = sqrt(hux)*sqrt(hxy)*sign(bux)*sign(bxy) * sqrt(huy)*sign(buy)
y2 = x2*bxy + scale(buy*u)[,1]*sqrt(huy) + scale(rnorm(n2, 0, 1))[,1]*sqrt(1-hxy-huy-2*cov_xy); # y = x + u + e;
y2 = scale(y2)[,1]
indexbuf = which((y2 >= qnorm(p, lower.tail=F)))
y2[indexbuf] = 1; y2[-indexbuf] = 0  # y, 0/1

# estimation
# x, simple regression, dataset 1
var_z1 = var(z1); var_z2 = var(z2)
bzx1_hat = cov(z1, x1)/var_z1; bzx1_hat_se = sqrt((1-var_z1*bzx1_hat^2)/(n1*var_z1))
bzx1_hat_pval = pchisq((bzx1_hat/bzx1_hat_se)^2, 1, lower.tail=F)
# y, logistic regression, dataset 2, WITHOUT selection bias
reg = summary(glm(y2~z2, family=binomial(link = "logit")))$coeff
bzy2_hat = reg[2,1]; bzy2_hat_se = reg[2,2]; bzy2_hat_pval = reg[2,4]; 

# dataset2, GWAS with selection
s_prob = gamma_0 + gamma_u*u + gamma_x*x2 # logit(pai) = gamma_0 + gamma_x + gamma_u + e
s_prob = exp(s_prob)/(1+exp(s_prob))
s = rbinom(n2, 1, s_prob) # selection event, 0/1
l1 = which(y2==0 & s==0); l2 = which(y2==1 & s==0)
z2_slct = z2[c(l1,l2)]; 
n2_slct = length(indexbuf)
y2_slct = y2[c(l1,l2)]
# y, logistic regression, dataset 2, WITH selection bias
reg = summary(glm(y2_slct~z2_slct, family=binomial(link = "logit")))$coeff
bzy2_slct_hat = reg[2,1]; bzy2_slct_hat_se = reg[2,2]; bzy2_slct_hat_pval = reg[2,4];

# SMR, WITHOUT selection bias
resbuf = SMR(bzx1_hat, bzx1_hat_se, bzx1_hat_pval, bzy2_hat, bzy2_hat_se)
bxy1_hat = resbuf$bxy; bxy1_hat_se = resbuf$bxy_se; bxy1_hat_pval = resbuf$bxy_pval
# SMR, WITH selection bias
resbuf = SMR(bzx1_hat, bzx1_hat_se, bzx1_hat_pval, bzy2_slct_hat, bzy2_slct_hat_se)
bxy2_hat = resbuf$bxy; bxy2_hat_se = resbuf$bxy_se; bxy2_hat_pval = resbuf$bxy_pval

# sample sizes
n1_case = length(which(y2==1)); n1_cntl = length(which(y2==0))
n2_case = length(l2); n2_cntl = length(l1)

# save the results
rst = rbind(cbind("Withou_slct", bxy1_hat, bxy1_hat_se, bxy1_hat_pval, n1, n1_case, n1_cntl, n1_case+n1_cntl),
            cbind("With_sclt", bxy2_hat, bxy2_hat_se, bxy2_hat_pval, n1, n2_case, n2_cntl, n2_case+n2_cntl))
colnames(rst) = c("slct", "bxy", "bxy_se", "bxy_pval", "n_x", "n_y_cases", "n_y_cntls", "n_y")
write.table(rst, rstfile, col.names=T, row.names=F, quote=F)
