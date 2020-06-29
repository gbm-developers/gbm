library("DALEX")

# titanic -----------------------------------------------------------------

titanic_data <- DALEX::titanic_imputed

titanic_glm_model <- glm(survived == 1 ~ class + gender + age + sibsp + 
                           parch + fare + embarked,
                         titanic_data, family = "binomial")

titanic_new_observation <- data.frame(
  class = factor("1st", levels = c("1st", "2nd", "3rd", "deck crew",
                                   "engineering crew", "restaurant staff",
                                   "victualling crew")),
  gender = factor("male", levels = c("female", "male")),
  age = 8,
  sibsp = 0,
  parch = 0,
  fare = 72,
  embarked = factor("Southampton", levels = c("Belfast","Cherbourg",
                                              "Queenstown","Southampton"))
)

titanic_aspects <- list(wealth = c("class", "fare"),
                        family = c("gender", "sibsp", "parch"),
                        age = "age",
                        embarked = "embarked")

# apartments --------------------------------------------------------------

apartments_data <- DALEX::apartments

apartments_lm_model <- lm(m2.price ~ ., data = apartments_data)
apartments_aspects <- list(space = c("surface", "no.rooms"),
                           construction.year = "construction.year",
                           floor = "floor",
                           district = "district")
apartments_new_observation <- DALEX::apartments_test[30,]
apartments_num <- apartments_data[,unlist(lapply(apartments_data, is.numeric))]
apartments_num_lm_model <- lm(m2.price ~ ., data = apartments_num)
apartments_num_new_observation <- apartments_num[30, ]
apartments_explain <- explain(model = apartments_num_lm_model,
                              data = apartments_num[,-1],
                              y = apartments_num[, 1],
                              verbose = FALSE)
