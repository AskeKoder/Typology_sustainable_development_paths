#Correlation based approach---------------------------------------
# #Try chats attempt at pred reduction
# reduce_pred_matrix <- function(df, pred, top_n = 10) {
#   num_vars <- names(df)[sapply(df, is.numeric)]
#   cors <- cor(df[num_vars], use = "pairwise.complete.obs")
#   
#   for (v in num_vars) {
#     if (v == "Prim_School_EnrollyearID2020"){next}
#     # exclude target itself
#     rel <- abs(cors[, v])
#     rel[v] <- 0
#     
#     # order predictors by absolute correlation
#     best <- names(sort(rel, decreasing = TRUE))[1:top_n]
#     
#     # reset row
#     pred[v, ] <- 0
#     pred[v, best] <- 1
#   }
#   return(pred)
# }
# 
# # Apply
# pred_reduced <- reduce_pred_matrix(data_wide, pred, top_n = 25)
# table(rowSums(pred_reduced))