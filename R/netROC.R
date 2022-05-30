#' @title Plot the combined ROC
#' @description This is a function to plot ROC of combined study.
#'
#' @param case_dir string.The directory of diseased data set or a single diseased file.
#' @param control_dir string.The directory of healthy data set or a single healthy file.
#' @param marker the combined markers identified by NetMoss.
#' @param metadata the metadata of all input studies.
#' @param plot.roc logical.If TURE then the combined ROC will be plotted.
#' @param train.num numerical. training times for the module. 20 was set by default.
#'
#' @return a TPR and FPR table of the combined study.
#' @export
#' @importFrom grDevices dev.off pdf
#'
#' @examples
#' data(testData)
#' case_dir = mydata[[1]]
#' control_dir = mydata[[2]]
#' nodes_result = NetMoss(case_dir = case_dir,
#'      control_dir = control_dir,
#'      net_case_dir = mydata[[3]],
#'      net_control_dir = mydata[[4]])
#' result = nodes_result[[1]]
#' marker = data.frame(result[which(result$p.adj < 0.05),])
#' marker = data.frame(marker[which(marker$NetMoss_Score > 0.3),])
#' rownames(marker) = marker$taxon_names
#' metadata = mydata[[5]]
#' myROC = netROC(case_dir = case_dir,
#'      control_dir = control_dir,
#'      marker = marker,
#'      metadata = metadata,
#'      plot.roc = TRUE,
#'      train.num = 20)
#'
netROC <-
  function(case_dir,
           control_dir,
           marker,
           metadata,
           plot.roc = TRUE,
           train.num = 20)
  {
    my.wd2 = getwd()

    message ("training datasets...")

    if (length(case_dir) > 1 || length(control_dir) > 1)
    {
      ####read
      x.case.net = case_dir
      if (colnames(x.case.net)[1] == "X")
      {
        rownames(x.case.net) = x.case.net[, 1]
      }

      x.control.net = control_dir
      if (colnames(x.control.net)[1] == "X")
      {
        rownames(x.control.net) = x.control.net[, 1]
      }
      rownames(metadata) = metadata$sample_id
      metadata$type = as.factor(metadata$type)

      ##loop
      roc.all = list()
      h = x.control.net[, -1]
      d = x.case.net[, -1]

      ##trim samples
      hh = h
      hh$genus = rownames(hh)
      hh = melt(hh, id.vars = "genus")
      dd = d
      dd$genus = rownames(dd)
      dd = melt(dd, id.vars = "genus")
      all1 = rbind(hh, dd)
      all2 = dcast(all1, genus ~ variable, mean, fill = 0)
      sample.all = all2[,-1]
      rownames(sample.all) = all2$genus
      aa = data.frame(t(sample.all))
      aa$sum = rowSums(aa)
      aa = aa[which(aa$sum != 0),]
      sample.all2 = data.frame(aa[,-ncol(aa)])
      m.marker = intersect(as.character(colnames(sample.all2)), as.character(rownames(marker)))
      m.meta = metadata[as.character(rownames(metadata)) %in% as.character(rownames(sample.all2)),]
      auc_crc_all = data.frame()
      my.marker = list()

      for (mm in 1:train.num)
      {
        ##train and test data
        train.data = sample(rownames(sample.all2), nrow(sample.all2) * 0.7, replace = F)
        sample.train = sample.all2[as.character(rownames(sample.all2)) %in% as.character(train.data),]
        sample.train$group = m.meta[as.character(rownames(sample.train)), "type"]
        test.data = setdiff(rownames(sample.all2), train.data)
        sample.test = sample.all2[as.character(rownames(sample.all2)) %in% as.character(test.data),]
        sample.test$group = m.meta[as.character(rownames(sample.test)), "type"]

        ##training
        if (length(intersect(m.marker, colnames(sample.train))) == 1)
        {
          sample.train2 = data.frame(sample.train[, as.character(colnames(sample.train)) %in% as.character(m.marker)])
          colnames(sample.train2) = as.character(intersect(m.marker, colnames(sample.train)))
          rownames(sample.train2) = rownames(sample.train)
        } else
        {
          sample.train2 = sample.train[, as.character(colnames(sample.train)) %in% as.character(m.marker)]
        }
        sample.train2$group = sample.train$group
        rf.train = randomForest(
          group ~ .,
          sample.train2,
          ntree = 1000,
          nperm = 100,
          importance = T
        )
        importance_rf = data.frame(importance(rf.train))
        importance_rf = importance_rf[order(importance_rf$MeanDecreaseAccuracy, decreasing = T),]

        #########10-fold cross validation
        train.5_10 = replicate(5,
                               rfcv(
                                 sample.train2[-ncol(sample.train2)],
                                 sample.train2$group,
                                 cv.fold = 10,
                                 step = 1.5
                               ),
                               simplify = F)
        train.5_10.2 = data.frame(sapply(train.5_10, '[[', 'error.cv'))
        train.5_10.2$names = rownames(train.5_10.2)
        train.5_10.2 = melt(train.5_10.2, id = 'names')
        train.5_10.2$names = as.numeric(as.character(train.5_10.2$names))
        train.5_10.2 = summaryBy(value ~ names, train.5_10.2, FUN = mean)

        ####re-trainning
        marker.num = min(train.5_10.2[which(train.5_10.2$value.mean == min(train.5_10.2$value.mean)), 1])
        marker.re = data.frame(rownames(importance_rf[1:marker.num,]))
        colnames(marker.re) = "Name"
        if (length(intersect(marker.re$Name, colnames(sample.train))) == 1)
        {
          sample.train3 = data.frame(sample.train[, as.character(colnames(sample.train)) %in% as.character(marker.re$Name)])
          colnames(sample.train3) = as.character(intersect(marker.re$Name, colnames(sample.train)))
          rownames(sample.train3) = rownames(sample.train)
        } else
        {
          sample.train3 = sample.train[, as.character(colnames(sample.train)) %in% as.character(marker.re$Name)]
        }
        sample.train3$group = sample.train$group
        rf.train = randomForest(
          group ~ .,
          sample.train3,
          ntree = 1000,
          nperm = 100,
          importance = T
        )

        ##marker
        marker.re2 = marker.re
        marker.re2$NetMoss_Score = marker[as.character(marker.re$Name),"NetMoss_Score"]
        my.marker[[mm]] = marker.re2

        ##testing
        test_pred.select = predict(rf.train, sample.test)
        test_freq.select = table(test_pred.select, sample.test$group)

        ##ROC
        group_select = data.frame(test_pred.select)
        group_select$test_orig.select = sample.test[as.character(rownames(group_select)), 'group']
        e11 = group_select
        e11$test_pred.select = as.numeric(e11$test_pred.select)
        e11$test_orig.select = as.numeric(e11$test_orig.select)
        e11[which(e11$test_orig.select != 1), 2] = 0
        e11[which(e11$test_pred.select != 1), 1] = 0
        for (i in 1:nrow(e11))
          ifelse(e11[i, 1] == e11[i, 2], e11[i, 3] <-
                   1, e11[i, 3] <- 0)
        b = data.frame(predict(rf.train, sample.test, type = 'prob'))
        e11$prob = b[as.character(rownames(e11)), "disease"]

        ##plot
        p = ggplot(e11, aes(d = test_orig.select, m = prob)) +
          geom_roc(n.cuts = 0) +
          style_roc()
        auc <- calc_auc(p)

        #roc
        eval(parse(text = paste("roc_crc", mm, " = e11", sep = "")))
        eval(parse(text = paste("auc_crc_all[", mm, ",1] = ", mm, sep = "")))
        eval(parse(text = paste(
          "auc_crc_all[", mm, ",2] = auc$AUC", sep = ""
        )))

        #pick
        if (auc$AUC > 0.7)
        {
          eval(parse(text = paste("roc.all[[1]] = roc_crc", mm, sep = "")))
          eval(parse(text = paste("my.marker.new = my.marker[[", mm, "]]",sep = "")))
          break
        }
        if (mm == train.num)
        {
          m.auc = which(auc_crc_all$V2 == max(auc_crc_all$V2))
          eval(parse(text = paste(
            "roc.all[[1]] = roc_crc", m.auc, sep = ""
          )))
          eval(parse(text = paste(
            "my.marker.new = my.marker[[", mm, "]]",sep = ""
          )))
        }
      }


      ################################
      mydata = roc.all[[1]]
      y <- mydata$test_orig.select
      y.hat <- mydata$prob
      tp.fp <- WeightedROC(y.hat, y)

      ##plot combined ROC
      setwd(my.wd2)

      if (plot.roc)
      {
        p2 =  ggplot() +
          geom_path(aes(FPR, TPR), data = tp.fp) +
          coord_equal() +
          annotate(
            "text",
            x = .75,
            y = .15,
            label = paste("AUC =", round(WeightedAUC(tp.fp), 2))
          ) +
          labs(x = "False positive fraction", y = "True positive fraction") +
          theme_bw()

        ggsave("NetMoss_ROC.pdf", p2)
        print (paste0("the ROC curve saved in ", my.wd2))
      }

      my.marker.new = my.marker.new[order(my.marker.new$NetMoss_Score,decreasing = T),]

      tp.fp2 = list()
      tp.fp2[[1]] = tp.fp[, c(3, 1, 2)]
      tp.fp2[[2]] = my.marker.new

      return(tp.fp2)

    } else if (file_test('-f', case_dir) &&
               file_test('-f', control_dir))
    {
      ####read
      x.case.net = read.table(
        file = case_dir,
        header = T,
        sep = "\t",
        row.names = 1
      )
      x.control.net = read.table(
        file = control_dir,
        header = T,
        sep = "\t",
        row.names = 1
      )

      ##loop
      roc.all = list()
      h = x.control.net
      d = x.case.net

      ##trim samples
      hh = h
      hh$genus = rownames(hh)
      hh = melt(hh, id.vars = "genus")
      dd = d
      dd$genus = rownames(dd)
      dd = melt(dd, id.vars = "genus")
      all1 = rbind(hh, dd)
      all2 = dcast(all1, genus ~ variable, mean, fill = 0)
      sample.all = all2[,-1]
      rownames(sample.all) = all2$genus
      aa = data.frame(t(sample.all))
      aa$sum = rowSums(aa)
      aa = aa[which(aa$sum != 0),]
      sample.all2 = data.frame(aa[,-ncol(aa)])
      m.marker = intersect(as.character(colnames(sample.all2)), as.character(rownames(marker)))
      m.meta = metadata[as.character(rownames(metadata)) %in% as.character(rownames(sample.all2)),]
      auc_crc_all = data.frame()
      my.marker = list()

      for (mm in 1:train.num)
      {
        ##train and test data
        train.data = sample(rownames(sample.all2), nrow(sample.all2) * 0.7, replace = F)
        sample.train = sample.all2[as.character(rownames(sample.all2)) %in% as.character(train.data),]
        sample.train$group = m.meta[as.character(rownames(sample.train)), "type"]
        test.data = setdiff(rownames(sample.all2), train.data)
        sample.test = sample.all2[as.character(rownames(sample.all2)) %in% as.character(test.data),]
        sample.test$group = m.meta[as.character(rownames(sample.test)), "type"]

        ##training
        if (length(intersect(m.marker, colnames(sample.train))) == 1)
        {
          sample.train2 = data.frame(sample.train[, as.character(colnames(sample.train)) %in% as.character(m.marker)])
          colnames(sample.train2) = as.character(intersect(m.marker, colnames(sample.train)))
          rownames(sample.train2) = rownames(sample.train)
        } else
        {
          sample.train2 = sample.train[, as.character(colnames(sample.train)) %in% as.character(m.marker)]
        }
        sample.train2$group = sample.train$group
        rf.train = randomForest(
          group ~ .,
          sample.train2,
          ntree = 1000,
          nperm = 100,
          importance = T
        )
        importance_rf = data.frame(importance(rf.train))
        importance_rf = importance_rf[order(importance_rf$MeanDecreaseAccuracy, decreasing = T),]

        #########10-fold cross validation
        train.5_10 = replicate(5,
                               rfcv(
                                 sample.train2[-ncol(sample.train2)],
                                 sample.train2$group,
                                 cv.fold = 10,
                                 step = 1.5
                               ),
                               simplify = F)
        train.5_10.2 = data.frame(sapply(train.5_10, '[[', 'error.cv'))
        train.5_10.2$names = rownames(train.5_10.2)
        train.5_10.2 = melt(train.5_10.2, id = 'names')
        train.5_10.2$names = as.numeric(as.character(train.5_10.2$names))
        train.5_10.2 = summaryBy(value ~ names, train.5_10.2, FUN = mean)

        ####re-trainning
        marker.num = min(train.5_10.2[which(train.5_10.2$value.mean == min(train.5_10.2$value.mean)), 1])
        marker.re = data.frame(rownames(importance_rf[1:marker.num,]))
        colnames(marker.re) = "Name"
        if (length(intersect(marker.re$Name, colnames(sample.train))) == 1)
        {
          sample.train3 = data.frame(sample.train[, as.character(colnames(sample.train)) %in% as.character(marker.re$Name)])
          colnames(sample.train3) = as.character(intersect(marker.re$Name, colnames(sample.train)))
          rownames(sample.train3) = rownames(sample.train)
        } else
        {
          sample.train3 = sample.train[, as.character(colnames(sample.train)) %in% as.character(marker.re$Name)]
        }
        sample.train3$group = sample.train$group
        rf.train = randomForest(
          group ~ .,
          sample.train3,
          ntree = 1000,
          nperm = 100,
          importance = T
        )

        ##marker
        marker.re2 = marker.re
        marker.re2$NetMoss_Score = marker[as.character(marker.re$Name),"NetMoss_Score"]
        my.marker[[mm]] = marker.re2

        ##testing
        test_pred.select = predict(rf.train, sample.test)
        test_freq.select = table(test_pred.select, sample.test$group)

        ##ROC
        group_select = data.frame(test_pred.select)
        group_select$test_orig.select = sample.test[as.character(rownames(group_select)), 'group']
        e11 = group_select
        e11$test_pred.select = as.numeric(e11$test_pred.select)
        e11$test_orig.select = as.numeric(e11$test_orig.select)
        e11[which(e11$test_orig.select != 1), 2] = 0
        e11[which(e11$test_pred.select != 1), 1] = 0
        for (i in 1:nrow(e11))
          ifelse(e11[i, 1] == e11[i, 2], e11[i, 3] <-
                   1, e11[i, 3] <- 0)
        b = data.frame(predict(rf.train, sample.test, type = 'prob'))
        e11$prob = b[as.character(rownames(e11)), "disease"]

        ##plot
        p = ggplot(e11, aes(d = test_orig.select, m = prob)) +
          geom_roc(n.cuts = 0) +
          style_roc()
        auc <- calc_auc(p)

        #roc
        eval(parse(text = paste("roc_crc", mm, " = e11", sep = "")))
        eval(parse(text = paste("auc_crc_all[", mm, ",1] = ", mm, sep = "")))
        eval(parse(text = paste(
          "auc_crc_all[", mm, ",2] = auc$AUC", sep = ""
        )))

        #pick
        if (auc$AUC > 0.7)
        {
          eval(parse(text = paste("roc.all[[1]] = roc_crc", mm, sep = "")))
          eval(parse(text = paste("my.marker.new = my.marker[[", mm, "]]",sep = "")))
          break
        }
        if (mm == train.num)
        {
          m.auc = which(auc_crc_all$V2 == max(auc_crc_all$V2))
          eval(parse(text = paste(
            "roc.all[[1]] = roc_crc", m.auc, sep = ""
          )))
          eval(parse(text = paste(
            "my.marker.new = my.marker[[", mm, "]]",sep = ""
          )))
        }

      }


      ################################
      mydata = roc.all[[1]]
      y <- mydata$test_orig.select
      y.hat <- mydata$prob
      tp.fp <- WeightedROC(y.hat, y)

      ##plot combined ROC
      setwd(my.wd2)

      if (plot.roc)
      {
        p2 =  ggplot() +
          geom_path(aes(FPR, TPR), data = tp.fp) +
          coord_equal() +
          annotate(
            "text",
            x = .75,
            y = .15,
            label = paste("AUC =", round(WeightedAUC(tp.fp), 2))
          ) +
          labs(x = "False positive fraction", y = "True positive fraction") +
          theme_bw()

        ggsave("NetMoss_ROC.pdf", p2)
        print (paste0("the ROC curve saved in ", my.wd2))
      }

      my.marker.new = my.marker.new[order(my.marker.new$NetMoss_Score,decreasing = T),]

      tp.fp2 = list()
      tp.fp2[[1]] = tp.fp[, c(3, 1, 2)]
      tp.fp2[[2]] = my.marker.new

      return(tp.fp2)

    } else
    {
      #case dir
      case_data_list <- list()
      case_samples <- c()
      setwd(case_dir)
      dir3 = list.files()
      n3 = length(dir3)
      for (m in 1:n3)
      {
        x.case.net = read.table(
          file = dir3[m],
          header = T,
          sep = "\t",
          row.names = 1
        )
        case_data_list[[m]] = x.case.net
        case_samples <-
          c(case_samples, length(colnames(x.case.net)) - 1)
      }
      #control dir
      control_data_list <- list()
      control_samples <- c()
      setwd(control_dir)
      dir4 = list.files()
      n4 = length(dir4)
      for (m in 1:n4)
      {
        x.control.net = read.table(
          file = dir4[m],
          header = T,
          sep = "\t",
          row.names = 1
        )
        control_data_list[[m]] = x.control.net
        control_samples <-
          c(control_samples, length(colnames(x.control.net)) - 1)
      }

      ##loop
      roc.all = list()
      my.marker.new = list()

      for (kk in 1:n3)
      {
        eval(parse(text = paste(
          "h = control_data_list[[", kk, "]]", sep = ""
        )))
        eval(parse(text = paste(
          "d = case_data_list[[", kk, "]]", sep = ""
        )))

        ##trim samples
        hh = h
        hh$genus = rownames(hh)
        hh = melt(hh, id.vars = "genus")
        dd = d
        dd$genus = rownames(dd)
        dd = melt(dd, id.vars = "genus")
        all1 = rbind(hh, dd)
        all2 = dcast(all1, genus ~ variable, mean, fill = 0)
        sample.all = all2[,-1]
        rownames(sample.all) = all2$genus
        aa = data.frame(t(sample.all))
        aa$sum = rowSums(aa)
        aa = aa[which(aa$sum != 0),]
        sample.all2 = data.frame(aa[,-ncol(aa)])
        m.marker = intersect(as.character(colnames(sample.all2)), as.character(rownames(marker)))
        m.meta = metadata[as.character(rownames(metadata)) %in% as.character(rownames(sample.all2)),]

        eval(parse(text = paste(
          "auc_crc", kk, "_all = data.frame()", sep = ""
        )))

        eval(parse(text = paste(
          "my.marker", kk, "_all = list()", sep = ""
        )))

        for (mm in 1:train.num)
        {
          ##train and test data
          train.data = sample(rownames(sample.all2),
                              nrow(sample.all2) * 0.7,
                              replace = F)
          sample.train = sample.all2[as.character(rownames(sample.all2)) %in% as.character(train.data),]
          sample.train$group = m.meta[as.character(rownames(sample.train)), "type"]
          test.data = setdiff(rownames(sample.all2), train.data)
          sample.test = sample.all2[as.character(rownames(sample.all2)) %in% as.character(test.data),]
          sample.test$group = m.meta[as.character(rownames(sample.test)), "type"]

          ##training
          if (length(intersect(m.marker, colnames(sample.train))) == 1)
          {
            sample.train2 = data.frame(sample.train[, as.character(colnames(sample.train)) %in% as.character(m.marker)])
            colnames(sample.train2) = as.character(intersect(m.marker, colnames(sample.train)))
            rownames(sample.train2) = rownames(sample.train)
          } else
          {
            sample.train2 = sample.train[, as.character(colnames(sample.train)) %in% as.character(m.marker)]
          }
          sample.train2$group = sample.train$group
          rf.train = randomForest(
            group ~ .,
            sample.train2,
            ntree = 1000,
            nperm = 100,
            importance = T
          )
          importance_rf = data.frame(importance(rf.train))
          importance_rf = importance_rf[order(importance_rf$MeanDecreaseAccuracy, decreasing = T),]

          #########10-fold cross validation
          train.5_10 = replicate(
            5,
            rfcv(
              sample.train2[-ncol(sample.train2)],
              sample.train2$group,
              cv.fold = 10,
              step = 1.5
            ),
            simplify = F
          )
          train.5_10.2 = data.frame(sapply(train.5_10, '[[', 'error.cv'))
          train.5_10.2$names = rownames(train.5_10.2)
          train.5_10.2 = melt(train.5_10.2, id = 'names')
          train.5_10.2$names = as.numeric(as.character(train.5_10.2$names))
          train.5_10.2 = summaryBy(value ~ names, train.5_10.2, FUN = mean)

          ####re-trainning
          marker.num = min(train.5_10.2[which(train.5_10.2$value.mean == min(train.5_10.2$value.mean)), 1])
          marker.re = data.frame(rownames(importance_rf[1:marker.num,]))
          colnames(marker.re) = "Name"
          if (length(intersect(marker.re$Name, colnames(sample.train))) == 1)
          {
            sample.train3 = data.frame(sample.train[, as.character(colnames(sample.train)) %in% as.character(marker.re$Name)])
            colnames(sample.train3) = as.character(intersect(marker.re$Name, colnames(sample.train)))
            rownames(sample.train3) = rownames(sample.train)
          } else
          {
            sample.train3 = sample.train[, as.character(colnames(sample.train)) %in% as.character(marker.re$Name)]
          }
          sample.train3$group = sample.train$group
          rf.train = randomForest(
            group ~ .,
            sample.train3,
            ntree = 1000,
            nperm = 100,
            importance = T
          )

          ##marker
          marker.re2 = marker.re
          marker.re2$NetMoss_Score = marker[as.character(marker.re$Name),"NetMoss_Score"]
          eval(parse(text = paste(
            "my.marker", kk, "_all[[",mm,"]] = marker.re2", sep = ""
          )))

          ##testing
          test_pred.select = predict(rf.train, sample.test)
          test_freq.select = table(test_pred.select, sample.test$group)

          ##ROC
          group_select = data.frame(test_pred.select)
          group_select$test_orig.select = sample.test[as.character(rownames(group_select)), 'group']
          e11 = group_select
          e11$test_pred.select = as.numeric(e11$test_pred.select)
          e11$test_orig.select = as.numeric(e11$test_orig.select)
          e11[which(e11$test_orig.select != 1), 2] = 0
          e11[which(e11$test_pred.select != 1), 1] = 0
          for (i in 1:nrow(e11))
            ifelse(e11[i, 1] == e11[i, 2], e11[i, 3] <-
                     1, e11[i, 3] <- 0)
          b = data.frame(predict(rf.train, sample.test, type = 'prob'))
          e11$prob = b[as.character(rownames(e11)), "disease"]

          ##plot
          p = ggplot(e11, aes(d = test_orig.select, m = prob)) +
            geom_roc(n.cuts = 0) +
            style_roc()
          auc <- calc_auc(p)

          #roc
          eval(parse(text = paste(
            "roc_crc", kk, "_", mm, " = e11", sep = ""
          )))
          eval(parse(text = paste(
            "auc_crc", kk, "_all[", mm, ",1] = ", mm, sep = ""
          )))
          eval(parse(
            text = paste("auc_crc", kk, "_all[", mm, ",2] = auc$AUC", sep = "")
          ))

          #pick
          if (auc$AUC > 0.7)
          {
            eval(parse(
              text = paste("roc.all[[", kk, "]] = roc_crc", kk, "_", mm, sep = "")
            ))
            eval(parse(
              text = paste("my.marker.new[[", kk, "]] = my.marker", kk, "_all[[",mm,"]]", sep = "")
            ))
            break
          }
          if (mm == train.num)
          {
            eval(parse(
              text = paste(
                "m.auc = which(auc_crc",
                kk,
                "_all$V2 == max(auc_crc",
                kk,
                "_all$V2))",
                sep = ""
              )
            ))
            eval(parse(
              text = paste("roc.all[[", kk, "]] = roc_crc", kk, "_", m.auc, sep = "")
            ))
            eval(parse(
              text = paste("my.marker.new[[", kk, "]] = my.marker", kk, "_all[[",mm,"]]", sep = "")
            ))
          }
        }

      }

      ##
      for (nn in 1:n3)
      {
        eval(parse(text = paste(
          "crc", nn, " = roc.all[[", nn, "]]", sep = ""
        )))
        eval(parse(text = paste(
          "crc", nn, "$type = 'CRC", nn, "'", sep = ""
        )))
      }

      if (n3 > 1)
      {
        mydata = crc1
        for (nn in 2:n3)
        {
          eval(parse(text = paste(
            "mydata = rbind(mydata,crc", nn, ")", sep = ""
          )))
        }

        ##weighted
        w.case = data.frame(sample = case_samples)
        w.case$w1 = 1 / w.case$sample * sum(w.case$sample)
        w.case$w2 = w.case$w1 / sum(w.case$w1)
        y <- mydata$test_orig.select
        w = rep(w.case[1, 3], nrow(crc1))
        for (nn in 2:n3)
        {
          eval(parse(
            text = paste(
              "w = c(w,rep(w.case[",
              nn,
              ",3],nrow(crc",
              nn,
              ")))",
              sep = ""
            )
          ))
        }
      } else
      {
        mydata = crc1

        ##weighted
        w.case = data.frame(sample = case_samples)
        w.case$w1 = 1 / w.case$sample * sum(w.case$sample)
        w.case$w2 = w.case$w1 / sum(w.case$w1)
        y <- mydata$test_orig.select
        w = rep(w.case[1, 3], nrow(crc1))
      }

      y.hat <- mydata$prob
      tp.fp <- WeightedROC(y.hat, y, w)

      ##plot combined ROC
      setwd(my.wd2)

      if (plot.roc)
      {
        p2 =  ggplot() +
          geom_path(aes(FPR, TPR), data = tp.fp) +
          coord_equal() +
          annotate(
            "text",
            x = .75,
            y = .15,
            label = paste("AUC =", round(WeightedAUC(tp.fp), 2))
          ) +
          labs(x = "False positive fraction", y = "True positive fraction") +
          theme_bw()

        ggsave("NetMoss_ROC.pdf", p2)
        print (paste0("the ROC curve saved in ", my.wd2))
      }

      ##marker
      all.marker = c()
      for (i in 1:n3)
      {
        all.marker = c(all.marker,as.character(my.marker.new[[i]]$Name))
      }
      all.marker2 = unique(all.marker)
      all.marker3 = data.frame(Name = all.marker2,
                               NetMoss_Score = marker[as.character(all.marker2),"NetMoss_Score"])
      all.marker3 = all.marker3[order(all.marker3$NetMoss_Score,decreasing = T),]

      tp.fp2 = list()
      tp.fp2[[1]] = tp.fp[, c(3, 1, 2)]
      tp.fp2[[2]] = all.marker3

      return(tp.fp2)

    }

  }
