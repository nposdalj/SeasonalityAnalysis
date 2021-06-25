#GAM w/ GEE

#GLMs

#GLM in days
hist(dayBinTAB$HoursProp)
glm_day = glm(HoursProp ~ day, data = dayBinTAB, family = poisson())
plot(glm_day)
summary(glm_day)

#other code
#fitting a simple GEE to view the output
gee_simp = GEE(HoursProp ~ day + Season, family = "poisson", data = oneyear)
summary(gee_simp)

#looks for statistical significane between seasons using a GLM
GLM_form = formula(dayBinTAB$HoursProp ~ dayBinTAB$Season)
gee1 = geeglm(GLM_form, data = dayBinTAB, id = Season, family=poisson("identity"))
gee1
coef(gee1)
vcov(gee1)
summary(gee1)


#GAM w/ GEE

#GLMs

#GLM in days
hist(dayBinTAB$HoursProp)
glm_day = glm(HoursProp ~ day, data = dayBinTAB, family = poisson())
plot(glm_day)
summary(glm_day)

#other code
#fitting a simple GEE to view the output
gee_simp = GEE(HoursProp ~ day + Season, family = "poisson", data = oneyear)
summary(gee_simp)

#looks for statistical significane between seasons using a GLM
GLM_form = formula(dayBinTAB$HoursProp ~ dayBinTAB$Season)
gee1 = geeglm(GLM_form, data = dayBinTAB, id = Season, family=poisson("identity"))
gee1
coef(gee1)
vcov(gee1)
summary(gee1)

