test_that("NetMoss() identifies markers", {
  data(testData)
  nodes_result = NetMoss(case_dir = mydata[[1]],
                         control_dir = mydata[[2]],
                         net_case_dir = mydata[[3]],
                         net_control_dir = mydata[[4]])
})
