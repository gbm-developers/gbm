pkgname <- "creditmodel"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('creditmodel')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("PCA_reduce")
### * PCA_reduce

flush(stderr()); flush(stdout())

### Name: PCA_reduce
### Title: PCA Dimension Reduction
### Aliases: PCA_reduce

### ** Examples

## Not run: 
##D num_x_list = get_names(dat = UCICreditCard, types = c('numeric'),
##D ex_cols = "ID$|date$|default.payment.next.month$", get_ex = FALSE)
##D  PCA_dat = PCA_reduce(train = UCICreditCard[num_x_list])
## End(Not run)



cleanEx()
nameEx("as_percent")
### * as_percent

flush(stderr()); flush(stdout())

### Name: as_percent
### Title: Percent Format
### Aliases: as_percent

### ** Examples

as_percent(0.2363, digits = 2)
as_percent(1)



cleanEx()
nameEx("char_cor_vars")
### * char_cor_vars

flush(stderr()); flush(stdout())

### Name: char_cor_vars
### Title: Cramer's V matrix between categorical variables.
### Aliases: char_cor_vars char_cor

### ** Examples

## Not run: 
##D char_x_list = get_names(dat = UCICreditCard,
##D types = c('factor', 'character'),
##D ex_cols = "ID$|date$|default.payment.next.month$", get_ex = FALSE)
##D  char_cor(dat = UCICreditCard[char_x_list])
## End(Not run)



cleanEx()
nameEx("char_to_num")
### * char_to_num

flush(stderr()); flush(stdout())

### Name: char_to_num
### Title: character to number
### Aliases: char_to_num

### ** Examples

dat_sub = lendingclub[c('dti_joint',	'emp_length')]
str(dat_sub)
#variables that are converted to numbers containing strings
dat_sub = char_to_num(dat_sub)
str(dat_sub)



cleanEx()
nameEx("check_rules")
### * check_rules

flush(stderr()); flush(stdout())

### Name: check_rules
### Title: check rules
### Aliases: check_rules

### ** Examples

train_test = train_test_split(UCICreditCard, split_type = "Random", prop = 0.8, save_data = FALSE)
dat_train = train_test$train
dat_test = train_test$test
dat_train$default.payment.next.month = as.numeric(dat_train$default.payment.next.month)
rules_list = get_ctree_rules(tree_fit = NULL, train_dat = dat_train[, 8:26],
                             target ="default.payment.next.month", test_dat = dat_test)[1:3,2]
check_rules(rules_list = rules_list, target = "default.payment.next.month",
test_dat = dat_test, value = NULL)




cleanEx()
nameEx("checking_data")
### * checking_data

flush(stderr()); flush(stdout())

### Name: checking_data
### Title: Checking Data
### Aliases: checking_data

### ** Examples

dat = checking_data(dat = UCICreditCard, target = "default.payment.next.month")



cleanEx()
nameEx("cor_heat_plot")
### * cor_heat_plot

flush(stderr()); flush(stdout())

### Name: cor_heat_plot
### Title: Correlation Heat Plot
### Aliases: cor_heat_plot

### ** Examples

train_test <- train_test_split(UCICreditCard,
split_type = "Random", prop = 0.8,save_data = FALSE)
dat_train = train_test$train
dat_test = train_test$test
cor_mat = cor(dat_train[,8:12],use = "complete.obs")
cor_heat_plot(cor_mat)



cleanEx()
nameEx("cor_plot")
### * cor_plot

flush(stderr()); flush(stdout())

### Name: cor_plot
### Title: Correlation Plot
### Aliases: cor_plot

### ** Examples

train_test <- train_test_split(UCICreditCard,
split_type = "Random", prop = 0.8,save_data = FALSE)
dat_train = train_test$train
dat_test = train_test$test
cor_plot(dat_train[,8:12],plot_show = TRUE)



cleanEx()
nameEx("cross_table")
### * cross_table

flush(stderr()); flush(stdout())

### Name: cross_table
### Title: cross_table
### Aliases: cross_table

### ** Examples

cross_table(dat = UCICreditCard, cross_x = "SEX",cross_y = "AGE",
 target = "default.payment.next.month", cross_type = "bad_pct",value = "LIMIT_BAL")
cross_table(dat = UCICreditCard, cross_x = c("SEX", "MARRIAGE"), cross_y = "AGE",
target = "default.payment.next.month", cross_type = "bad_pct",value = "LIMIT_BAL")



cleanEx()
nameEx("customer_segmentation")
### * customer_segmentation

flush(stderr()); flush(stdout())

### Name: customer_segmentation
### Title: Customer Segmentation
### Aliases: customer_segmentation

### ** Examples

clust <- customer_segmentation(dat = lendingclub[1:10000,20:30],
                              x_list = NULL, ex_cols = "id$|loan_status",
                              cluster_control = list(meth = "FCM", kc = 2),  save_data = FALSE,
                              tree_control = list(minbucket = round(nrow(lendingclub) / 10)),
                              file_name = NULL, dir_path = tempdir())



cleanEx()
nameEx("cut_equal")
### * cut_equal

flush(stderr()); flush(stdout())

### Name: cut_equal
### Title: Generating Initial Equal Size Sample Bins
### Aliases: cut_equal

### ** Examples

#equal sample size breaks
equ_breaks = cut_equal(dat = UCICreditCard[, "PAY_AMT2"], g = 10)




cleanEx()
nameEx("cv_split")
### * cv_split

flush(stderr()); flush(stdout())

### Name: cv_split
### Title: Stratified Folds
### Aliases: cv_split

### ** Examples

sub = cv_split(UCICreditCard, k = 30)[[1]]
dat = UCICreditCard[sub,]



cleanEx()
nameEx("data_cleansing")
### * data_cleansing

flush(stderr()); flush(stdout())

### Name: data_cleansing
### Title: Data Cleaning
### Aliases: data_cleansing

### ** Examples

#data cleaning
dat_cl = data_cleansing(dat = UCICreditCard[1:2000,],
                       target = "default.payment.next.month",
                       x_list = NULL,
                       obs_id = "ID",
                       occur_time = "apply_date",
                       ex_cols = c("PAY_6|BILL_"),
                       outlier_proc = TRUE,
                       missing_proc = TRUE,
                       low_var = TRUE,
                       save_data = FALSE)




cleanEx()
nameEx("data_exploration")
### * data_exploration

flush(stderr()); flush(stdout())

### Name: data_exploration
### Title: Data Exploration
### Aliases: data_exploration

### ** Examples

data_ex = data_exploration(dat = UCICreditCard[1:1000,])



cleanEx()
nameEx("date_cut")
### * date_cut

flush(stderr()); flush(stdout())

### Name: date_cut
### Title: Date Time Cut Point
### Aliases: date_cut

### ** Examples

date_cut(dat_time = lendingclub$issue_d, pct = 0.8)
#"2018-08-01"



cleanEx()
nameEx("de_one_hot_encoding")
### * de_one_hot_encoding

flush(stderr()); flush(stdout())

### Name: de_one_hot_encoding
### Title: Recovery One-Hot Encoding
### Aliases: de_one_hot_encoding

### ** Examples

#one hot encoding
dat1 = one_hot_encoding(dat = UCICreditCard,
cat_vars = c("SEX", "MARRIAGE"),
merge_cat = TRUE, na_act = TRUE)
#de one hot encoding
dat2 = de_one_hot_encoding(dat_one_hot = dat1,
cat_vars = c("SEX","MARRIAGE"),
na_act = FALSE)



cleanEx()
nameEx("de_percent")
### * de_percent

flush(stderr()); flush(stdout())

### Name: de_percent
### Title: Recovery Percent Format
### Aliases: de_percent

### ** Examples

de_percent("24%")



cleanEx()
nameEx("digits_num")
### * digits_num

flush(stderr()); flush(stdout())

### Name: digits_num
### Title: Number of digits
### Aliases: digits_num

### ** Examples

## Not run: 
##D digits_num(lendingclub[,"dti"])
##D # 7
## End(Not run)



cleanEx()
nameEx("entropy_weight")
### * entropy_weight

flush(stderr()); flush(stdout())

### Name: entropy_weight
### Title: Entropy Weight Method
### Aliases: entropy_weight

### ** Examples

entropy_weight(dat = ewm_data,ID = "ID",
              pos_vars = -c(7,11),
              neg_vars = c(7,11))



cleanEx()
nameEx("entry_rate_na")
### * entry_rate_na

flush(stderr()); flush(stdout())

### Name: entry_rate_na
### Title: Max Percent of missing Value
### Aliases: entry_rate_na

### ** Examples

datss = entry_rate_na(dat = lendingclub[1:1000, ], nr = 0.98)



cleanEx()
nameEx("fast_high_cor_filter")
### * fast_high_cor_filter

flush(stderr()); flush(stdout())

### Name: fast_high_cor_filter
### Title: high_cor_filter
### Aliases: fast_high_cor_filter high_cor_filter

### ** Examples

# calculate iv for each variable.
iv_list = feature_selector(dat_train = UCICreditCard[1:1000,], dat_test = NULL,
target = "default.payment.next.month",
occur_time = "apply_date",
filter = c("IV"), cv_folds = 1, iv_cp = 0.01,
ex_cols = "ID$|date$|default.payment.next.month$",
save_data = FALSE, vars_name = FALSE)
fast_high_cor_filter(dat = UCICreditCard[1:1000,],
com_list = iv_list, save_data = FALSE,
ex_cols = "ID$|date$|default.payment.next.month$",
p = 0.9, cor_class = FALSE ,var_name = FALSE)



cleanEx()
nameEx("feature_selector")
### * feature_selector

flush(stderr()); flush(stdout())

### Name: feature_selector
### Title: Feature Selection Wrapper
### Aliases: feature_selector

### ** Examples

feature_selector(dat_train = UCICreditCard[1:1000,c(2,8:12,26)],
                      dat_test = NULL, target = "default.payment.next.month",
                      occur_time = "apply_date", filter = c("IV", "PSI"),
                      cv_folds = 1, iv_cp = 0.01, psi_cp = 0.1, xgb_cp = 0, cor_cp = 0.98,
                      vars_name = FALSE,note = FALSE)



cleanEx()
nameEx("gather_data")
### * gather_data

flush(stderr()); flush(stdout())

### Name: gather_data
### Title: gather or aggregate data
### Aliases: gather_data

### ** Examples

dat = data.frame(id = c(1,1,1,2,2,3,3,3,4,4,4,4,4,5,5,6,7,7,
                            8,8,8,9,9,9,10,10,11,11,11,11,11,11),
                     terms = c('a','b','c','a','c','d','d','a',
                               'b','c','a','c','d','a','c',
                                  'd','a','e','f','b','c','f','b',
                               'c','h','h','i','c','d','g','k','k'),
                     time = c(8,3,1,9,6,1,4,9,1,3,4,8,2,7,1,
                              3,4,1,8,7,2,5,7,8,8,2,1,5,7,2,7,3))

gather_data(dat = dat, x_list = "time", ID = 'id', FUN = sum_x)



cleanEx()
nameEx("gbm_filter")
### * gbm_filter

flush(stderr()); flush(stdout())

### Name: gbm_filter
### Title: Select Features using GBM
### Aliases: gbm_filter

### ** Examples

GBM.params = gbm_params(n.trees = 2, interaction.depth = 2, shrinkage = 0.1,
                       bag.fraction = 1, train.fraction = 1,
                       n.minobsinnode = 30,
                     cv.folds = 2)
## Not run: 
##D  features <- gbm_filter(dat = UCICreditCard[1:1000, c(8:12, 26)],
##D          target = "default.payment.next.month",
##D       occur_time = "apply_date",
##D      GBM.params = GBM.params
##D        , vars_name = FALSE)
## End(Not run)



cleanEx()
nameEx("get_bins_table_all")
### * get_bins_table_all

flush(stderr()); flush(stdout())

### Name: get_bins_table_all
### Title: Table of Binning
### Aliases: get_bins_table_all get_bins_table

### ** Examples

breaks_list = get_breaks_all(dat = UCICreditCard, x_list = names(UCICreditCard)[3:4],
target = "default.payment.next.month", equal_bins =TRUE,best = FALSE,g=5,
ex_cols = "ID|apply_date", save_data = FALSE)
get_bins_table_all(dat = UCICreditCard, breaks_list = breaks_list,
target = "default.payment.next.month")



cleanEx()
nameEx("get_breaks_all")
### * get_breaks_all

flush(stderr()); flush(stdout())

### Name: get_breaks_all
### Title: Generates Best Breaks for Binning
### Aliases: get_breaks_all get_breaks

### ** Examples

#controls
tree_control = list(p = 0.02, cp = 0.000001, xval = 5, maxdepth = 10)
bins_control = list(bins_num = 10, bins_pct = 0.02, b_chi = 0.02, b_odds = 0.1,
                   b_psi = 0.05, b_or = 15, mono = 0.2, odds_psi = 0.1, kc = 5)
# get categrory variable breaks
b <-  get_breaks(dat = UCICreditCard[1:1000,], x = "MARRIAGE",
                target = "default.payment.next.month",
                occur_time = "apply_date",
                sp_values = list(-1, "missing"),
                tree_control = tree_control, bins_control = bins_control)
# get numeric variable breaks
b2 <-  get_breaks(dat = UCICreditCard[1:1000,], x = "PAY_2",
                 target = "default.payment.next.month",
                 occur_time = "apply_date",
                 sp_values = list(-1, "missing"),
                 tree_control = tree_control, bins_control = bins_control)
# get breaks of all predictive variables
b3 <-  get_breaks_all(dat = UCICreditCard[1:1000,], target = "default.payment.next.month",
                     x_list = c("MARRIAGE","PAY_2"),
                     occur_time = "apply_date", ex_cols = "ID",
                     sp_values = list(-1, "missing"),
                    tree_control = tree_control, bins_control = bins_control,
                     save_data = FALSE)




cleanEx()
nameEx("get_correlation_group")
### * get_correlation_group

flush(stderr()); flush(stdout())

### Name: get_correlation_group
### Title: get_correlation_group
### Aliases: get_correlation_group select_cor_group select_cor_list

### ** Examples

## Not run: 
##D cor_mat = cor(UCICreditCard[8:20],
##D use = "complete.obs", method = "spearman")
##D get_correlation_group(cor_mat, p = 0.6 )
## End(Not run)



cleanEx()
nameEx("get_ctree_rules")
### * get_ctree_rules

flush(stderr()); flush(stdout())

### Name: get_ctree_rules
### Title: Parse desision tree rules
### Aliases: get_ctree_rules

### ** Examples

train_test <- train_test_split(UCICreditCard, split_type = "Random", prop = 0.8, save_data = FALSE)
dat_train = train_test$train
dat_test = train_test$test
dat_train$default.payment.next.month = as.numeric(dat_train$default.payment.next.month)
get_ctree_rules(tree_fit = NULL, train_dat = dat_train[, 8:26],
target ="default.payment.next.month", test_dat = dat_test)




cleanEx()
nameEx("get_iv_all")
### * get_iv_all

flush(stderr()); flush(stdout())

### Name: get_iv_all
### Title: Calculate Information Value (IV) 'get_iv' is used to calculate
###   Information Value (IV) of an independent variable. 'get_iv_all' can
###   loop through IV for all specified independent variables.
### Aliases: get_iv_all get_iv

### ** Examples

get_iv_all(dat = UCICreditCard,
 x_list = names(UCICreditCard)[3:10],
 equal_bins = TRUE, best = FALSE,
 target = "default.payment.next.month",
 ex_cols = "ID|apply_date")
get_iv(UCICreditCard, x = "PAY_3",
       equal_bins = TRUE, best = FALSE,
 target = "default.payment.next.month")



cleanEx()
nameEx("get_logistic_coef")
### * get_logistic_coef

flush(stderr()); flush(stdout())

### Name: get_logistic_coef
### Title: get logistic coef
### Aliases: get_logistic_coef

### ** Examples

# dataset spliting
sub = cv_split(UCICreditCard, k = 30)[[1]]
dat = UCICreditCard[sub,]
#rename the target variable
dat = re_name(dat, "default.payment.next.month", "target")
dat = data_cleansing(dat, target = "target", obs_id = "ID",
occur_time = "apply_date", miss_values =  list("", -1))
#train_ test pliting
train_test <- train_test_split(dat, split_type = "OOT", prop = 0.7,
                                occur_time = "apply_date")
dat_train = train_test$train
dat_test = train_test$test
#get breaks of all predictive variables
x_list = c("PAY_0", "LIMIT_BAL", "PAY_AMT5", "EDUCATION", "PAY_3", "PAY_2")
breaks_list <- get_breaks_all(dat = dat_train, target = "target",
                              x_list = x_list, occur_time = "apply_date", ex_cols = "ID",
save_data = FALSE, note = FALSE)
#woe transforming
train_woe = woe_trans_all(dat = dat_train,
                          target = "target",
                          breaks_list = breaks_list,
                          woe_name = FALSE)
test_woe = woe_trans_all(dat = dat_test,
                       target = "target",
                         breaks_list = breaks_list,
                         note = FALSE)
Formula = as.formula(paste("target", paste(x_list, collapse = ' + '), sep = ' ~ '))
set.seed(46)
lr_model = glm(Formula, data = train_woe[, c("target", x_list)], family = binomial(logit))
#get LR coefficient
dt_imp_LR = get_logistic_coef(lg_model = lr_model, save_data = FALSE)
bins_table = get_bins_table_all(dat = dat_train, target = "target",
                                x_list = x_list,dat_test = dat_test,
                               breaks_list = breaks_list, note = FALSE)
#score card
LR_score_card <- get_score_card(lg_model = lr_model, bins_table, target = "target")
#scoring
train_pred = dat_train[, c("ID", "apply_date", "target")]
test_pred = dat_test[, c("ID", "apply_date", "target")]
train_pred$pred_LR = score_transfer(model = lr_model,
                                                    tbl_woe = train_woe,
                                                    save_data = TRUE)[, "score"]

test_pred$pred_LR = score_transfer(model = lr_model,
tbl_woe = test_woe, save_data = FALSE)[, "score"]



cleanEx()
nameEx("get_names")
### * get_names

flush(stderr()); flush(stdout())

### Name: get_names
### Title: Get Variable Names
### Aliases: get_names

### ** Examples

x_list = get_names(dat = UCICreditCard, types = c('factor', 'character'),
ex_cols = c("default.payment.next.month","ID$|_date$"), get_ex = FALSE)
x_list = get_names(dat = UCICreditCard, types = c('numeric', 'character', "integer"),
ex_cols = c("default.payment.next.month", "ID$|SEX "), get_ex = FALSE)



cleanEx()
nameEx("get_plots")
### * get_plots

flush(stderr()); flush(stdout())

### Name: get_plots
### Title: Plot Independent Variables Distribution
### Aliases: get_plots plot_vars

### ** Examples

train_test <- train_test_split(UCICreditCard[1:1000,], split_type = "Random",
 prop = 0.8, save_data = FALSE)
dat_train = train_test$train
dat_test = train_test$test
get_plots(dat_train[, c(8, 26)], dat_test = dat_test[, c(8, 26)],
target = "default.payment.next.month")



cleanEx()
nameEx("get_psi_all")
### * get_psi_all

flush(stderr()); flush(stdout())

### Name: get_psi_all
### Title: Calculate Population Stability Index (PSI) 'get_psi' is used to
###   calculate Population Stability Index (PSI) of an independent
###   variable. 'get_psi_all' can loop through PSI for all specified
###   independent variables.
### Aliases: get_psi_all get_psi

### ** Examples

#  dat_test is null
get_psi(dat = UCICreditCard, x = "PAY_3", occur_time = "apply_date")
# dat_test is not all
# train_test split
train_test = train_test_split(dat = UCICreditCard, prop = 0.7, split_type = "OOT",
                             occur_time = "apply_date", start_date = NULL, cut_date = NULL,
                            save_data = FALSE, note = FALSE)
dat_ex = train_test$train
dat_ac = train_test$test
# generate psi table
get_psi(dat = dat_ex, dat_test = dat_ac, x = "PAY_3",
       occur_time = "apply_date", bins_no = TRUE)



cleanEx()
nameEx("get_psi_iv_all")
### * get_psi_iv_all

flush(stderr()); flush(stdout())

### Name: get_psi_iv_all
### Title: Calculate IV & PSI
### Aliases: get_psi_iv_all get_psi_iv

### ** Examples

iv_list = get_psi_iv_all(dat = UCICreditCard[1:1000, ],
x_list = names(UCICreditCard)[3:5], equal_bins = TRUE,
target = "default.payment.next.month", ex_cols = "ID|apply_date")
get_psi_iv(UCICreditCard, x = "PAY_3",
target = "default.payment.next.month",bins_total = TRUE)



cleanEx()
nameEx("get_psi_plots")
### * get_psi_plots

flush(stderr()); flush(stdout())

### Name: get_psi_plots
### Title: Plot PSI(Population Stability Index)
### Aliases: get_psi_plots psi_plot

### ** Examples

train_test <- train_test_split(UCICreditCard[1:1000,], split_type = "Random",
 prop = 0.8, save_data = FALSE)
dat_train = train_test$train
dat_test = train_test$test
get_psi_plots(dat_train[, c(8, 9)], dat_test = dat_test[, c(8, 9)])



cleanEx()
nameEx("get_score_card")
### * get_score_card

flush(stderr()); flush(stdout())

### Name: get_score_card
### Title: Score Card
### Aliases: get_score_card

### ** Examples

# dataset spliting
sub = cv_split(UCICreditCard, k = 30)[[1]]
dat = UCICreditCard[sub,]
#rename the target variable
dat = re_name(dat, "default.payment.next.month", "target")
dat = data_cleansing(dat, target = "target", obs_id = "ID",
occur_time = "apply_date", miss_values =  list("", -1))
#train_ test pliting
train_test <- train_test_split(dat, split_type = "OOT", prop = 0.7,
                                occur_time = "apply_date")
dat_train = train_test$train
dat_test = train_test$test
#get breaks of all predictive variables
x_list = c("PAY_0", "LIMIT_BAL", "PAY_AMT5", "EDUCATION", "PAY_3", "PAY_2")
breaks_list <- get_breaks_all(dat = dat_train, target = "target",
                              x_list = x_list, occur_time = "apply_date", ex_cols = "ID",
save_data = FALSE, note = FALSE)
#woe transforming
train_woe = woe_trans_all(dat = dat_train,
                          target = "target",
                          breaks_list = breaks_list,
                          woe_name = FALSE)
test_woe = woe_trans_all(dat = dat_test,
                       target = "target",
                         breaks_list = breaks_list,
                         note = FALSE)
Formula = as.formula(paste("target", paste(x_list, collapse = ' + '), sep = ' ~ '))
set.seed(46)
lr_model = glm(Formula, data = train_woe[, c("target", x_list)], family = binomial(logit))
#get LR coefficient
dt_imp_LR = get_logistic_coef(lg_model = lr_model, save_data = FALSE)
bins_table = get_bins_table_all(dat = dat_train, target = "target",
                                 dat_test = dat_test,
                                x_list = x_list,
                               breaks_list = breaks_list, note = FALSE)
#score card
LR_score_card <- get_score_card(lg_model = lr_model, bins_table, target = "target")
#scoring
train_pred = dat_train[, c("ID", "apply_date", "target")]
test_pred = dat_test[, c("ID", "apply_date", "target")]
train_pred$pred_LR = score_transfer(model = lr_model,
                                                    tbl_woe = train_woe,
                                                    save_data = FALSE)[, "score"]

test_pred$pred_LR = score_transfer(model = lr_model,
tbl_woe = test_woe, save_data = FALSE)[, "score"]



cleanEx()
nameEx("get_tree_breaks")
### * get_tree_breaks

flush(stderr()); flush(stdout())

### Name: get_tree_breaks
### Title: Getting the breaks for terminal nodes from decision tree
### Aliases: get_tree_breaks

### ** Examples

#tree breaks
tree_control = list(p = 0.02, cp = 0.000001, xval = 5, maxdepth = 10)
tree_breaks = get_tree_breaks(dat = UCICreditCard, x = "MARRIAGE",
target = "default.payment.next.month", tree_control = tree_control)



cleanEx()
nameEx("get_x_list")
### * get_x_list

flush(stderr()); flush(stdout())

### Name: get_x_list
### Title: Get X List.
### Aliases: get_x_list

### ** Examples

x_list = get_x_list(x_list = NULL,dat_train = UCICreditCard,
ex_cols = c("default.payment.next.month","ID$|_date$"))



cleanEx()
nameEx("grapes-alike-grapes")
### * grapes-alike-grapes

flush(stderr()); flush(stdout())

### Name: %alike%
### Title: Fuzzy String matching
### Aliases: %alike%

### ** Examples

"xyz"  %alike% "xy"



cleanEx()
nameEx("grapes-islike-grapes")
### * grapes-islike-grapes

flush(stderr()); flush(stdout())

### Name: %islike%
### Title: Fuzzy String matching
### Aliases: %islike%

### ** Examples

 "xyz"  %islike% "yz$"



cleanEx()
nameEx("is_date")
### * is_date

flush(stderr()); flush(stdout())

### Name: is_date
### Title: is_date
### Aliases: is_date

### ** Examples

is_date(lendingclub$issue_d)



cleanEx()
nameEx("ks_table")
### * ks_table

flush(stderr()); flush(stdout())

### Name: ks_table
### Title: ks_table & plot
### Aliases: ks_table ks_table_plot ks_psi_plot model_key_index

### ** Examples

sub = cv_split(UCICreditCard, k = 30)[[1]]
dat = UCICreditCard[sub,]
dat = re_name(dat, "default.payment.next.month", "target")
dat = data_cleansing(dat, target = "target", obs_id = "ID",
occur_time = "apply_date", miss_values = list("", -1))

train_test <- train_test_split(dat, split_type = "OOT", prop = 0.7,
                                occur_time = "apply_date")
dat_train = train_test$train
dat_test = train_test$test
x_list = c("PAY_0", "LIMIT_BAL", "PAY_AMT5", "PAY_3", "PAY_2")
Formula = as.formula(paste("target", paste(x_list, collapse = ' + '), sep = ' ~ '))
set.seed(46)
lr_model = glm(Formula, data = dat_train[, c("target", x_list)], family = binomial(logit))

dat_train$pred_LR = round(predict(lr_model, dat_train[, x_list], type = "response"), 5)
dat_test$pred_LR = round(predict(lr_model, dat_test[, x_list], type = "response"), 5)
# model evaluation
ks_psi_plot(train_pred = dat_train, test_pred = dat_test,
                            score = "pred_LR", target = "target",
                            plot_show = TRUE)
tb_pred <- ks_table_plot(train_pred = dat_train, test_pred = dat_test,
                                        score = "pred_LR", target = "target",
                                     g = 10, g_width = 13, plot_show = FALSE)
key_index = model_key_index(tb_pred)



cleanEx()
nameEx("lasso_filter")
### * lasso_filter

flush(stderr()); flush(stdout())

### Name: lasso_filter
### Title: Variable selection by LASSO
### Aliases: lasso_filter

### ** Examples

 sub = cv_split(UCICreditCard, k = 40)[[1]]
 dat = UCICreditCard[sub,]
 dat = re_name(dat, "default.payment.next.month", "target")
 dat_train = data_cleansing(dat, target = "target", obs_id = "ID", occur_time = "apply_date",
  miss_values = list("", -1))
 dat_train = process_nas(dat_train)
 #get breaks of all predictive variables
 x_list = c("PAY_0", "LIMIT_BAL", "PAY_AMT5", "EDUCATION", "PAY_3", "PAY_2")
 breaks_list <- get_breaks_all(dat = dat_train, target = "target",
                                x_list = x_list, occur_time = "apply_date", ex_cols = "ID",
  save_data = FALSE, note = FALSE)
 #woe transform
 train_woe = woe_trans_all(dat = dat_train,
                            target = "target",
                            breaks_list = breaks_list,
                            woe_name = FALSE)
 lasso_filter(dat_train = train_woe, 
         target = "target", x_list = x_list,
       save_data = FALSE, plot.it = FALSE)



cleanEx()
nameEx("log_trans")
### * log_trans

flush(stderr()); flush(stdout())

### Name: log_trans
### Title: Logarithmic transformation
### Aliases: log_trans log_vars

### ** Examples

dat = log_trans(dat = UCICreditCard, target = "default.payment.next.month",
x_list =NULL,cor_dif = 0.01,ex_cols = "ID", note = TRUE)



cleanEx()
nameEx("loop_function")
### * loop_function

flush(stderr()); flush(stdout())

### Name: loop_function
### Title: Loop Function. #' 'loop_function' is an iterator to loop through
### Aliases: loop_function

### ** Examples

dat = UCICreditCard[24:26]
num_x_list = get_names(dat = dat, types = c('numeric', 'integer', 'double'),
                      ex_cols = NULL, get_ex = FALSE)
dat[ ,num_x_list] = loop_function(func = outliers_kmeans_lof, x_list = num_x_list,
                                   args = list(dat = dat),
                                   bind = "cbind", as_list = FALSE,
                                 parallel = FALSE)



cleanEx()
nameEx("love_color")
### * love_color

flush(stderr()); flush(stdout())

### Name: love_color
### Title: love_color
### Aliases: love_color

### ** Examples

love_color(color="dark_cyan")



cleanEx()
nameEx("low_variance_filter")
### * low_variance_filter

flush(stderr()); flush(stdout())

### Name: low_variance_filter
### Title: Filtering Low Variance Variables
### Aliases: low_variance_filter

### ** Examples

dat = low_variance_filter(lendingclub[1:1000, ], lvp = 0.9)




cleanEx()
nameEx("lr_vif")
### * lr_vif

flush(stderr()); flush(stdout())

### Name: lr_vif
### Title: Variance-Inflation Factors
### Aliases: lr_vif

### ** Examples

sub = cv_split(UCICreditCard, k = 30)[[1]]
x_list = c("PAY_0", "LIMIT_BAL", "PAY_AMT5", "PAY_3", "PAY_2")
dat = re_name(UCICreditCard[sub,], "default.payment.next.month", "target")
dat = dat[,c("target",x_list)]

dat = data_cleansing(dat, miss_values = list("", -1))

train_test <- train_test_split(dat,  prop = 0.7)
dat_train = train_test$train
dat_test = train_test$test

Formula = as.formula(paste("target", paste(x_list, collapse = ' + '), sep = ' ~ '))
set.seed(46)
lr_model = glm(Formula, data = dat_train[, c("target", x_list)], family = binomial(logit))
lr_vif(lr_model)
get_logistic_coef(lr_model)
class(dat)
mod = lr_model
lr_vif(lr_model)



cleanEx()
nameEx("max_min_norm")
### * max_min_norm

flush(stderr()); flush(stdout())

### Name: max_min_norm
### Title: Max Min Normalization
### Aliases: max_min_norm

### ** Examples

dat_s = apply(UCICreditCard[,12:14], 2, max_min_norm)



cleanEx()
nameEx("merge_category")
### * merge_category

flush(stderr()); flush(stdout())

### Name: merge_category
### Title: Merge Category
### Aliases: merge_category

### ** Examples

#merge_catagory
dat =  merge_category(lendingclub,ex_cols = "id$|_d$")
char_list = get_names(dat = dat,types = c('factor', 'character'),
ex_cols = "id$|_d$", get_ex = FALSE)
str(dat[,char_list])



cleanEx()
nameEx("min_max_norm")
### * min_max_norm

flush(stderr()); flush(stdout())

### Name: min_max_norm
### Title: Min Max Normalization
### Aliases: min_max_norm

### ** Examples

dat_s = apply(UCICreditCard[,12:14], 2, min_max_norm)



cleanEx()
nameEx("model_result_plot")
### * model_result_plot

flush(stderr()); flush(stdout())

### Name: model_result_plot
### Title: model result plots 'model_result_plot' is a wrapper of
###   following: 'perf_table' is for generating a model performance table.
###   'ks_plot' is for K-S. 'roc_plot' is for ROC. 'lift_plot' is for Lift
###   Chart. 'score_distribution_plot' is for ploting the score
###   distribution.
### Aliases: model_result_plot perf_table ks_plot lift_plot roc_plot
###   score_distribution_plot

### ** Examples

sub = cv_split(UCICreditCard, k = 30)[[1]]
dat = UCICreditCard[sub,]
dat = re_name(dat, "default.payment.next.month", "target")
x_list = c("PAY_0", "LIMIT_BAL", "PAY_AMT5", "PAY_3", "PAY_2")
dat = data_cleansing(dat, target = "target", obs_id = "ID",x_list = x_list,
occur_time = "apply_date", miss_values = list("", -1))
dat = process_nas(dat,default_miss = TRUE)
train_test <- train_test_split(dat, split_type = "OOT", prop = 0.7,
                                occur_time = "apply_date")
dat_train = train_test$train
dat_test = train_test$test
Formula = as.formula(paste("target", paste(x_list, collapse = ' + '), sep = ' ~ '))
set.seed(46)
lr_model = glm(Formula, data = dat_train[, c("target", x_list)], family = binomial(logit))

dat_train$pred_LR = round(predict(lr_model, dat_train[, x_list], type = "response"), 5)
dat_test$pred_LR = round(predict(lr_model, dat_test[, x_list], type = "response"), 5)
# model evaluation
perf_table(train_pred = dat_train, test_pred = dat_test, target = "target", score = "pred_LR")
ks_plot(train_pred = dat_train, test_pred = dat_test, target = "target", score = "pred_LR")
roc_plot(train_pred = dat_train, test_pred = dat_test, target = "target", score = "pred_LR")
#lift_plot(train_pred = dat_train, test_pred = dat_test, target = "target", score = "pred_LR")
#score_distribution_plot(train_pred = dat_train, test_pred = dat_test,
#target = "target", score = "pred_LR")
#model_result_plot(train_pred = dat_train, test_pred = dat_test,
#target = "target", score = "pred_LR")



cleanEx()
nameEx("multi_grid")
### * multi_grid

flush(stderr()); flush(stdout())

### Name: multi_grid
### Title: Arrange list of plots into a grid
### Aliases: multi_grid

### ** Examples

library(ggplot2)
sub = cv_split(UCICreditCard, k = 30)[[1]]
dat = UCICreditCard[sub,]
dat = re_name(dat, "default.payment.next.month", "target")
dat = data_cleansing(dat, target = "target", obs_id = "ID",
occur_time = "apply_date", miss_values = list("", -1))
dat = process_nas(dat)
train_test <- train_test_split(dat, split_type = "OOT", prop = 0.7,
                                occur_time = "apply_date")
dat_train = train_test$train
dat_test = train_test$test
x_list = c("PAY_0", "LIMIT_BAL", "PAY_AMT5", "PAY_3", "PAY_2")
Formula = as.formula(paste("target", paste(x_list, collapse = ' + '), sep = ' ~ '))
set.seed(46)
lr_model = glm(Formula, data = dat_train[, c("target", x_list)], family = binomial(logit))

dat_train$pred_LR = round(predict(lr_model, dat_train[, x_list], type = "response"), 5)
dat_test$pred_LR = round(predict(lr_model, dat_test[, x_list], type = "response"), 5)
# model evaluation
p1 =  ks_plot(train_pred = dat_train, test_pred = dat_test, target = "target", score = "pred_LR")
p2 =  roc_plot(train_pred = dat_train, test_pred = dat_test, target = "target", score = "pred_LR")
p3 =  lift_plot(train_pred = dat_train, test_pred = dat_test, target = "target", score = "pred_LR")
p4 = score_distribution_plot(train_pred = dat_train, test_pred = dat_test,
target = "target", score = "pred_LR")
p_plots= multi_grid(p1,p2,p3,p4)
plot(p_plots)



cleanEx()
nameEx("multi_left_join")
### * multi_left_join

flush(stderr()); flush(stdout())

### Name: multi_left_join
### Title: multi_left_join
### Aliases: multi_left_join

### ** Examples

multi_left_join(UCICreditCard[1:10, 1:10], UCICreditCard[1:10, c(1,8:14)],
UCICreditCard[1:10, c(1,20:25)], by = "ID")



cleanEx()
nameEx("null_blank_na")
### * null_blank_na

flush(stderr()); flush(stdout())

### Name: null_blank_na
### Title: Encode NAs
### Aliases: null_blank_na

### ** Examples

datss = null_blank_na(dat = UCICreditCard[1:1000, ], miss_values =list(-1,-2))



cleanEx()
nameEx("one_hot_encoding")
### * one_hot_encoding

flush(stderr()); flush(stdout())

### Name: one_hot_encoding
### Title: One-Hot Encoding
### Aliases: one_hot_encoding

### ** Examples

dat1 = one_hot_encoding(dat = UCICreditCard,
cat_vars = c("SEX", "MARRIAGE"),
merge_cat = TRUE, na_act = TRUE)
dat2 = de_one_hot_encoding(dat_one_hot = dat1,
cat_vars = c("SEX","MARRIAGE"), na_act = FALSE)




cleanEx()
nameEx("partial_dependence_plot")
### * partial_dependence_plot

flush(stderr()); flush(stdout())

### Name: partial_dependence_plot
### Title: partial_dependence_plot
### Aliases: partial_dependence_plot get_partial_dependence_plots

### ** Examples

sub = cv_split(UCICreditCard, k = 30)[[1]]
dat = UCICreditCard[sub,]
dat = re_name(dat, "default.payment.next.month", "target")
dat = data_cleansing(dat, target = "target", obs_id = "ID",
occur_time = "apply_date", miss_values = list("", -1))

train_test <- train_test_split(dat, split_type = "OOT", prop = 0.7,
                                occur_time = "apply_date")
dat_train = train_test$train
dat_test = train_test$test
x_list = c("PAY_0", "LIMIT_BAL", "PAY_AMT5", "PAY_3", "PAY_2")
Formula = as.formula(paste("target", paste(x_list, collapse = ' + '), sep = ' ~ '))
set.seed(46)
lr_model = glm(Formula, data = dat_train[, c("target", x_list)], family = binomial(logit))
#plot partial dependency of one variable
partial_dependence_plot(model = lr_model, x ="LIMIT_BAL", x_train = dat_train)
#plot partial dependency of all variables
pd_list = get_partial_dependence_plots(model = lr_model, x_list = x_list[1:2],
 x_train = dat_train, save_data = FALSE,plot_show = TRUE)



cleanEx()
nameEx("plot_bar")
### * plot_bar

flush(stderr()); flush(stdout())

### Name: plot_bar
### Title: Plot Bar
### Aliases: plot_bar

### ** Examples

plot_bar(dat = lendingclub, x = "grade")



cleanEx()
nameEx("plot_box")
### * plot_box

flush(stderr()); flush(stdout())

### Name: plot_box
### Title: Plot Box
### Aliases: plot_box

### ** Examples

plot_box(lendingclub, x = "grade", y = "installment", g = 7)



cleanEx()
nameEx("plot_density")
### * plot_density

flush(stderr()); flush(stdout())

### Name: plot_density
### Title: Plot Density
### Aliases: plot_density

### ** Examples

plot_density(dat = lendingclub, x = "annual_inc",y = "emp_length", m =0, hist = FALSE)
plot_density(dat = lendingclub, x = "annual_inc", m = 2,
colors_y = love_color(type = "line")[c(1,3)])



cleanEx()
nameEx("plot_distribution")
### * plot_distribution

flush(stderr()); flush(stdout())

### Name: plot_distribution
### Title: Plot Distribution
### Aliases: plot_distribution plot_distribution_x

### ** Examples

plot_distribution_x(dat = lendingclub, x = "max_bal_bc", g = 10,
	cut_bin = 'equal_width')
plot_distribution(dat = lendingclub, x_list = c("max_bal_bc", "installment"), 
     g = 10,dir_path = tempdir(),
	cut_bin = 'equal_width')



cleanEx()
nameEx("plot_oot_perf")
### * plot_oot_perf

flush(stderr()); flush(stdout())

### Name: plot_oot_perf
### Title: plot_oot_perf 'plot_oot_perf' is for ploting performance of
###   cross time samples in the future
### Aliases: plot_oot_perf

### ** Examples

sub = cv_split(UCICreditCard, k = 30)[[1]]
dat = UCICreditCard[sub,]
dat = re_name(dat, "default.payment.next.month", "target")
x_list = c("PAY_0", "LIMIT_BAL", "PAY_AMT5", "PAY_3", "PAY_2")
dat = data_cleansing(dat, target = "target", obs_id = "ID",x_list = x_list,
occur_time = "apply_date", miss_values = list("", -1))
dat = process_nas(dat)
train_test <- train_test_split(dat, split_type = "OOT", prop = 0.7,
                                occur_time = "apply_date")
dat_train = train_test$train
dat_test = train_test$test
Formula = as.formula(paste("target", paste(x_list, collapse = ' + '), sep = ' ~ '))
set.seed(46)
lr_model = glm(Formula, data = dat_train[, c("target", x_list)], family = binomial(logit))

dat_train$pred_LR = round(predict(lr_model, dat_train[, x_list], type = "response"), 5)
dat_test$pred_LR = round(predict(lr_model, dat_test[, x_list], type = "response"), 5)
plot_oot_perf(dat_test = dat_test, occur_time = "apply_date", target = "target", x = "pred_LR")



cleanEx()
nameEx("plot_relative_freq_histogram")
### * plot_relative_freq_histogram

flush(stderr()); flush(stdout())

### Name: plot_relative_freq_histogram
### Title: Plot Relative Frequency Histogram
### Aliases: plot_relative_freq_histogram

### ** Examples

plot_relative_freq_histogram(dat = lendingclub, x = "grade", y = "dti_joint", g = 7,
	cut_bin = 'equal_width')



cleanEx()
nameEx("plot_table")
### * plot_table

flush(stderr()); flush(stdout())

### Name: plot_table
### Title: plot_table
### Aliases: plot_table

### ** Examples

iv_list = get_psi_iv_all(dat = UCICreditCard[1:1000, ],
                         x_list = names(UCICreditCard)[3:5], equal_bins = TRUE,
                         target = "default.payment.next.month", ex_cols = "ID|apply_date")
iv_dt =get_psi_iv(UCICreditCard, x = "PAY_3",
                  target = "default.payment.next.month", bins_total = TRUE)

plot_table(iv_dt)



cleanEx()
nameEx("process_nas")
### * process_nas

flush(stderr()); flush(stdout())

### Name: process_nas
### Title: missing Treatment
### Aliases: process_nas process_nas_var

### ** Examples

dat_na = process_nas(dat = UCICreditCard[1:1000,],
parallel = FALSE,ex_cols = "ID$", method = "median")




cleanEx()
nameEx("process_outliers")
### * process_outliers

flush(stderr()); flush(stdout())

### Name: process_outliers
### Title: Outliers Treatment
### Aliases: process_outliers outliers_kmeans_lof

### ** Examples

dat_out = process_outliers(UCICreditCard[1:10000,c(18:21,26)],
                        target = "default.payment.next.month",
                       ex_cols = "date$", kc = 3, kn = 10, 
                       parallel = FALSE,note = TRUE)



cleanEx()
nameEx("psi_iv_filter")
### * psi_iv_filter

flush(stderr()); flush(stdout())

### Name: psi_iv_filter
### Title: Variable reduction based on Information Value & Population
###   Stability Index filter
### Aliases: psi_iv_filter

### ** Examples

psi_iv_filter(dat= UCICreditCard[1:1000,c(2,4,8:9,26)],
             target = "default.payment.next.month",
             occur_time = "apply_date",
             parallel = FALSE)



cleanEx()
nameEx("quick_as_df")
### * quick_as_df

flush(stderr()); flush(stdout())

### Name: quick_as_df
### Title: List as data.frame quickly
### Aliases: quick_as_df

### ** Examples


UCICreditCard = quick_as_df(UCICreditCard)




cleanEx()
nameEx("ranking_percent_proc")
### * ranking_percent_proc

flush(stderr()); flush(stdout())

### Name: ranking_percent_proc
### Title: Ranking Percent Process
### Aliases: ranking_percent_proc ranking_percent_proc_x
###   ranking_percent_dict ranking_percent_dict_x

### ** Examples

rank_dict = ranking_percent_dict(dat = UCICreditCard[1:1000,],
x_list = c("LIMIT_BAL","BILL_AMT2","PAY_AMT3"), ex_cols = NULL )
UCICreditCard_new = ranking_percent_proc(dat = UCICreditCard[1:1000,],
x_list = c("LIMIT_BAL", "BILL_AMT2", "PAY_AMT3"), rank_dict = rank_dict, parallel = FALSE)



cleanEx()
nameEx("re_code")
### * re_code

flush(stderr()); flush(stdout())

### Name: re_code
### Title: re_code 're_code' search for matches to argument pattern within
###   each element of a character vector:
### Aliases: re_code

### ** Examples

SEX  = sample(c("F","M"),1000,replace = TRUE)
codes= data.frame(ori_value = c('F','M'), code = c(0,1) )
SEX_re = re_code(SEX,codes)



cleanEx()
nameEx("re_name")
### * re_name

flush(stderr()); flush(stdout())

### Name: re_name
### Title: Rename
### Aliases: re_name

### ** Examples

dt = re_name(dat = UCICreditCard, "default.payment.next.month" , "target")
names(dt['target'])



cleanEx()
nameEx("remove_duplicated")
### * remove_duplicated

flush(stderr()); flush(stdout())

### Name: remove_duplicated
### Title: Remove Duplicated Observations
### Aliases: remove_duplicated

### ** Examples

datss = remove_duplicated(dat = UCICreditCard,
target = "default.payment.next.month",
obs_id = "ID", occur_time =  "apply_date")



cleanEx()
nameEx("require_packages")
### * require_packages

flush(stderr()); flush(stdout())

### Name: require_packages
### Title: Packages required and intallment
### Aliases: require_packages

### ** Examples

## Not run: 
##D require_packages(data.table, ggplot2, dplyr)
## End(Not run)



cleanEx()
nameEx("rowAny")
### * rowAny

flush(stderr()); flush(stdout())

### Name: rowAny
### Title: Functions for vector operation.
### Aliases: rowAny rowAllnas colAllnas colAllzeros rowAll rowCVs rowSds
###   colSds rowMaxs rowMins rowMaxMins colMaxMins cnt_x sum_x max_x min_x
###   avg_x

### ** Examples

#any row has missing values
row_amy =  rowAny(UCICreditCard[8:10])
#rows which is all missing values
row_na =  rowAllnas(UCICreditCard[8:10])
#cols which is all missing values
col_na =  colAllnas(UCICreditCard[8:10])
#cols which is all zeros
row_zero =  colAllzeros(UCICreditCard[8:10])
#sum all numbers of a row
row_all =  rowAll(UCICreditCard[8:10])
#caculate cv of a row
row_cv =  rowCVs(UCICreditCard[8:10])
#caculate sd of a row
row_sd =  rowSds(UCICreditCard[8:10])
#caculate sd of a column
col_sd =  colSds(UCICreditCard[8:10])



cleanEx()
nameEx("rules_filter")
### * rules_filter

flush(stderr()); flush(stdout())

### Name: rules_filter
### Title: rules_filter
### Aliases: rules_filter

### ** Examples

train_test <- train_test_split(UCICreditCard, split_type = "Random", prop = 0.8, save_data = FALSE)
dat_train = train_test$train
dat_test = train_test$test
dat_train$default.payment.next.month = as.numeric(dat_train$default.payment.next.month)
rules_list = get_ctree_rules(tree_fit = NULL, train_dat = dat_train[, 8:26],
                             target ="default.payment.next.month", test_dat = dat_test)[1:3,2]
 new_dat = rules_filter(rules_list = rules_list[3], dat = dat_test)




cleanEx()
nameEx("rules_result")
### * rules_result

flush(stderr()); flush(stdout())

### Name: rules_result
### Title: rules_result
### Aliases: rules_result

### ** Examples

train_test <- train_test_split(UCICreditCard, split_type = "Random", prop = 0.8, save_data = FALSE)
dat_train = train_test$train
dat_test = train_test$test
dat_train$default.payment.next.month = as.numeric(dat_train$default.payment.next.month)
rules_list = get_ctree_rules(tree_fit = NULL, train_dat = dat_train[, 8:26],
                             target ="default.payment.next.month", test_dat = dat_test)[1:3,2]
dat_test$rules_result = rules_result(rules_list = rules_list[3], dat = dat_test)




cleanEx()
nameEx("save_data")
### * save_data

flush(stderr()); flush(stdout())

### Name: save_data
### Title: Save data
### Aliases: save_data

### ** Examples

save_data(UCICreditCard,"UCICreditCard", tempdir())



cleanEx()
nameEx("score_transfer")
### * score_transfer

flush(stderr()); flush(stdout())

### Name: score_transfer
### Title: Score Transformation
### Aliases: score_transfer

### ** Examples

# dataset spliting
sub = cv_split(UCICreditCard, k = 30)[[1]]
dat = UCICreditCard[sub,]
#rename the target variable
dat = re_name(dat, "default.payment.next.month", "target")
dat = data_cleansing(dat, target = "target", obs_id = "ID",
occur_time = "apply_date", miss_values =  list("", -1))
#train_ test pliting
train_test <- train_test_split(dat, split_type = "OOT", prop = 0.7,
                                occur_time = "apply_date")
dat_train = train_test$train
dat_test = train_test$test
#get breaks of all predictive variables
x_list = c("PAY_0", "LIMIT_BAL", "PAY_AMT5", "EDUCATION", "PAY_3", "PAY_2")
breaks_list <- get_breaks_all(dat = dat_train, target = "target",
                              x_list = x_list, occur_time = "apply_date", ex_cols = "ID",
save_data = FALSE, note = FALSE)
#woe transforming
train_woe = woe_trans_all(dat = dat_train,
                          target = "target",
                          breaks_list = breaks_list,
                          woe_name = FALSE)
test_woe = woe_trans_all(dat = dat_test,
                       target = "target",
                         breaks_list = breaks_list,
                         note = FALSE)
Formula = as.formula(paste("target", paste(x_list, collapse = ' + '), sep = ' ~ '))
set.seed(46)
lr_model = glm(Formula, data = train_woe[, c("target", x_list)], family = binomial(logit))
#get LR coefficient
dt_imp_LR = get_logistic_coef(lg_model = lr_model, save_data = FALSE)
bins_table = get_bins_table_all(dat = dat_train, target = "target",
                                x_list = x_list,dat_test = dat_test,
                               breaks_list = breaks_list, note = FALSE)
#score card
LR_score_card <- get_score_card(lg_model = lr_model, bins_table, target = "target")
#scoring
train_pred = dat_train[, c("ID", "apply_date", "target")]
test_pred = dat_test[, c("ID", "apply_date", "target")]
train_pred$pred_LR = score_transfer(model = lr_model,
                                                    tbl_woe = train_woe,
                                                    save_data = FALSE)[, "score"]

test_pred$pred_LR = score_transfer(model = lr_model,
tbl_woe = test_woe, save_data = FALSE)[, "score"]



cleanEx()
nameEx("select_best_class")
### * select_best_class

flush(stderr()); flush(stdout())

### Name: select_best_class
### Title: Generates Best Binning Breaks
### Aliases: select_best_class select_best_breaks

### ** Examples

#equal sample size breaks
equ_breaks = cut_equal(dat = UCICreditCard[, "PAY_AMT2"], g = 10)

# select best bins
bins_control = list(bins_num = 10, bins_pct = 0.02, b_chi = 0.02,
b_odds = 0.1, b_psi = 0.05, b_or = 0.15, mono = 0.3, odds_psi = 0.1, kc = 1)
select_best_breaks(dat = UCICreditCard, x = "PAY_AMT2", breaks = equ_breaks,
target = "default.payment.next.month", occur_time = "apply_date",
sp_values = NULL, bins_control = bins_control)



cleanEx()
nameEx("split_bins")
### * split_bins

flush(stderr()); flush(stdout())

### Name: split_bins
### Title: split_bins
### Aliases: split_bins

### ** Examples

bins = split_bins(dat = UCICreditCard,
x = "PAY_AMT1", breaks = NULL, bins_no = TRUE)



cleanEx()
nameEx("split_bins_all")
### * split_bins_all

flush(stderr()); flush(stdout())

### Name: split_bins_all
### Title: Split bins all
### Aliases: split_bins_all

### ** Examples

sub = cv_split(UCICreditCard, k = 30)[[1]]
dat = UCICreditCard[sub,]
dat = re_name(dat, "default.payment.next.month", "target")
dat = data_cleansing(dat, target = "target", obs_id = "ID", occur_time = "apply_date",
miss_values =  list("", -1))

train_test <- train_test_split(dat, split_type = "OOT", prop = 0.7,
                                occur_time = "apply_date")
dat_train = train_test$train
dat_test = train_test$test
#get breaks of all predictive variables
x_list = c("PAY_0", "LIMIT_BAL", "PAY_AMT5", "EDUCATION", "PAY_3", "PAY_2")
breaks_list <- get_breaks_all(dat = dat_train, target = "target",
                              x_list = x_list, occur_time = "apply_date", ex_cols = "ID",
save_data = FALSE, note  = FALSE)
#woe transform
train_bins = split_bins_all(dat = dat_train,
                          breaks_list = breaks_list,
                          woe_name = FALSE)
test_bins = split_bins_all(dat = dat_test,
                         breaks_list = breaks_list,
                         note = FALSE)




cleanEx()
nameEx("str_match")
### * str_match

flush(stderr()); flush(stdout())

### Name: str_match
### Title: string match #' 'str_match' search for matches to argument
###   pattern within each element of a character vector:
### Aliases: str_match

### ** Examples

orignal_nam = c("12mdd","11mdd","10mdd")
str_match(str_r = orignal_nam,pattern= "\\d+")



cleanEx()
nameEx("swap_analysis")
### * swap_analysis

flush(stderr()); flush(stdout())

### Name: swap_analysis
### Title: Swap Out/Swap In Analysis
### Aliases: swap_analysis

### ** Examples

swap_analysis(dat = UCICreditCard, new_rules = list("SEX == 'male' & AGE < 25"),
 old_rules = list("SEX == 'male' & AGE < 30"),
 target = "default.payment.next.month", cross_type = "bad_pct", value = "LIMIT_BAL")



cleanEx()
nameEx("term_tfidf")
### * term_tfidf

flush(stderr()); flush(stdout())

### Name: term_tfidf
### Title: TF-IDF
### Aliases: term_tfidf term_idf term_filter

### ** Examples

term_df = data.frame(id = c(1,1,1,2,2,3,3,3,4,4,4,4,4,5,5,6,7,7,
                            8,8,8,9,9,9,10,10,11,11,11,11,11,11),
terms = c('a','b','c','a','c','d','d','a','b','c','a','c','d','a','c',
          'd','a','e','f','b','c','f','b','c','h','h','i','c','d','g','k','k'))
term_df = term_filter(term_df = term_df, low_freq = 1)
idf = term_idf(term_df)
tf_idf = term_tfidf(term_df,idf = idf)



cleanEx()
nameEx("time_series_proc")
### * time_series_proc

flush(stderr()); flush(stdout())

### Name: time_series_proc
### Title: Process time series data
### Aliases: time_series_proc

### ** Examples

dat = data.frame(id = c(1,1,1,2,2,3,3,3,4,4,4,4,4,5,5,6,7,7,
                            8,8,8,9,9,9,10,10,11,11,11,11,11,11),
                     terms = c('a','b','c','a','c','d','d','a',
                               'b','c','a','c','d','a','c',
                                  'd','a','e','f','b','c','f','b',
                               'c','h','h','i','c','d','g','k','k'),
                     time = c(8,3,1,9,6,1,4,9,1,3,4,8,2,7,1,
                              3,4,1,8,7,2,5,7,8,8,2,1,5,7,2,7,3))

time_series_proc(dat = dat, ID = 'id', group = 'terms',time = 'time')



cleanEx()
nameEx("time_transfer")
### * time_transfer

flush(stderr()); flush(stdout())

### Name: time_transfer
### Title: Time Format Transfering
### Aliases: time_transfer

### ** Examples

#transfer a variable.
dat = time_transfer(dat = lendingclub,date_cols = "issue_d")
class(dat[,"issue_d"])
#transfer a group of variables with similar name.
dat = time_transfer(dat = lendingclub,date_cols = "_d$")
class(dat[,"issue_d"])
#transfer all time variables.
dat = time_transfer(dat = lendingclub,date_cols = NULL)
class(dat[,"issue_d"])



cleanEx()
nameEx("train_test_split")
### * train_test_split

flush(stderr()); flush(stdout())

### Name: train_test_split
### Title: Train-Test-Split
### Aliases: train_test_split

### ** Examples

train_test <- train_test_split(lendingclub,
split_type = "OOT", prop = 0.7,
occur_time = "issue_d", seed = 12, save_data = FALSE)
dat_train = train_test$train
dat_test = train_test$test



cleanEx()
nameEx("training_model")
### * training_model

flush(stderr()); flush(stdout())

### Name: training_model
### Title: Training model
### Aliases: training_model

### ** Examples

sub = cv_split(UCICreditCard, k = 30)[[1]]
dat = UCICreditCard[sub,]
x_list = c("LIMIT_BAL")
B_model = training_model(dat = dat,
                         model_name = "UCICreditCard",
                         target = "default.payment.next.month",
							x_list = x_list,
                         occur_time =NULL,
                         obs_id =NULL,
							dat_test = NULL,
                         preproc = FALSE,
                         outlier_proc = FALSE,
                         missing_proc = FALSE,
                         feature_filter = NULL,
                         algorithm = list("LR"),
                         LR.params = lr_params(lasso = FALSE,
                                               step_wise = FALSE,
                                                 score_card = FALSE),
                         breaks_list = NULL,
                         parallel = FALSE,
                         cores_num = NULL,
                         save_pmml = FALSE,
                         plot_show = FALSE,
                         vars_plot = FALSE,
                         model_path = tempdir(),
                         seed = 46)




cleanEx()
nameEx("var_group_proc")
### * var_group_proc

flush(stderr()); flush(stdout())

### Name: var_group_proc
### Title: Process group numeric variables
### Aliases: var_group_proc

### ** Examples

dat = data.frame(id = c(1,1,1,2,2,3,3,3,4,4,4,4,4,5,5,6,7,7,
                            8,8,8,9,9,9,10,10,11,11,11,11,11,11),
                     terms = c('a','b','c','a','c','d','d','a',
                               'b','c','a','c','d','a','c',
                                  'd','a','e','f','b','c','f','b',
                               'c','h','h','i','c','d','g','k','k'),
                     time = c(8,3,1,9,6,1,4,9,1,3,4,8,2,7,1,
                              3,4,1,8,7,2,5,7,8,8,2,1,5,7,2,7,3))

time_series_proc(dat = dat, ID = 'id', group = 'terms',time = 'time')



cleanEx()
nameEx("woe_trans_all")
### * woe_trans_all

flush(stderr()); flush(stdout())

### Name: woe_trans_all
### Title: WOE Transformation
### Aliases: woe_trans_all woe_trans

### ** Examples

sub = cv_split(UCICreditCard, k = 30)[[1]]
dat = UCICreditCard[sub,]
dat = re_name(dat, "default.payment.next.month", "target")
dat = data_cleansing(dat, target = "target", obs_id = "ID", occur_time = "apply_date",
miss_values =  list("", -1))

train_test <- train_test_split(dat, split_type = "OOT", prop = 0.7,
                                occur_time = "apply_date")
dat_train = train_test$train
dat_test = train_test$test
#get breaks of all predictive variables
x_list = c("PAY_0", "LIMIT_BAL", "PAY_AMT5", "EDUCATION", "PAY_3", "PAY_2")
breaks_list <- get_breaks_all(dat = dat_train, target = "target",
                              x_list = x_list, occur_time = "apply_date", ex_cols = "ID",
save_data = FALSE, note  = FALSE)
#woe transform
train_woe = woe_trans_all(dat = dat_train,
                          target = "target",
                          breaks_list = breaks_list,
                          woe_name = FALSE)
test_woe = woe_trans_all(dat = dat_test,
                       target = "target",
                         breaks_list = breaks_list,
                         note = FALSE)




cleanEx()
nameEx("xgb_filter")
### * xgb_filter

flush(stderr()); flush(stdout())

### Name: xgb_filter
### Title: Select Features using XGB
### Aliases: xgb_filter

### ** Examples

dat = UCICreditCard[1:1000,c(2,4,8:9,26)]
xgb_params = list(nrounds = 100, max_depth = 6, eta = 0.1,
                                       min_child_weight = 1, subsample = 1,
                                       colsample_bytree = 1, gamma = 0, scale_pos_weight = 1,
                                       early_stopping_rounds = 10,
                                       objective = "binary:logistic")
## Not run: 
##D xgb_features <- xgb_filter(dat_train = dat, dat_test = NULL,
##D target = "default.payment.next.month", occur_time = "apply_date",f_eval = 'ks',
##D xgb_params = xgb_params,
##D cv_folds = 1, ex_cols = "ID$|date$|default.payment.next.month$", vars_name = FALSE)
## End(Not run)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
