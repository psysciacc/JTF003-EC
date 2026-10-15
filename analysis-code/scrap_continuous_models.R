
# fit models for H1
mH1_all_cgi <- glmer(correct_phase2 ~ (feedback * justification) * CGI + 
                   (feedback + feedback:justification||Country)+ (1|problem), 
                 family = binomial("logit"), data = dat_H1, 
                 control = glmerControl(optimizer = "bobyqa", 
                                        optCtrl = list(maxfun = 1e5)))

# 
m_all_cgi <- glmer(corrected ~ feedback * justification * CGI + 
                 (feedback * justification||Country)+ (1|problem), 
               family = binomial("logit"), 
               data = dat_Hother, 
               control = glmerControl(optimizer = "bobyqa", 
                                      optCtrl = list(maxfun = 1e5)))

summary(mH1_all_cgi)
summary(m_all_cgi)
