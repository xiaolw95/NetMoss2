#' @title Integrating Networks.
#' @description This is a function to integrate data set from different studies based on network.
#'
#' @param case_dir string.The directory of diseased data set or a single diseased file.
#' @param control_dir string.The directory of healthy data set or a single healthy file.
#' @param net_case_dir string.The directory of network correlation of diseased data set or a single diseased network file.
#' @param net_control_dir string.The directory of network correlation of healthy data set or a single healthy network file.
#'
#' @return integrated diseased and healthy networks.
#' @export
#' @importFrom utils file_test read.table write.table
#'
#' @examples
#' data(testData)
#' netAll = getNetwork(case_dir = mydata[[1]],
#'      control_dir = mydata[[2]],
#'      net_case_dir = mydata[[3]],
#'      net_control_dir = mydata[[4]])
#'
getNetwork <-
  function(case_dir,
           control_dir,
           net_case_dir,
           net_control_dir) {
    message ("importing datasets...")

    ######net_case_dir is a data.frame#######
    if (length(net_case_dir) > 1 || length(net_control_dir) > 1)
    {
      ###read
      case_union = net_case_dir
      control_union = net_control_dir

      ##union
      union_genus1 <- rownames(case_union)
      union_genus2 <- rownames(control_union)

      union_genus <- union(union_genus1, union_genus2)

      union_matrix <-
        matrix(nrow = length(union_genus),
               ncol = length(union_genus),
               0)

      rownames(union_matrix) <- union_genus
      colnames(union_matrix) <- union_genus
      diag(union_matrix) <- 1

      ##combine
      #case
      fi <- case_union
      re_union <- union_matrix
      diff_genus <- setdiff(union_genus, rownames(fi))
      rownames(re_union) <- c(rownames(fi), diff_genus)
      colnames(re_union) <- c(rownames(fi), diff_genus)
      for (j in 1:nrow(fi)) {
        re_union[j, 1:nrow(fi)] <- as.numeric(fi[j,])
      }
      re_union <- re_union[union_genus, union_genus]
      diag(re_union) <- 1
      case_union_data_list <- re_union

      #control
      fi <- control_union
      re_union <- union_matrix
      diff_genus <- setdiff(union_genus, rownames(fi))
      rownames(re_union) <- c(rownames(fi), diff_genus)
      colnames(re_union) <- c(rownames(fi), diff_genus)
      for (j in 1:nrow(fi)) {
        re_union[j, 1:nrow(fi)] <- as.numeric(fi[j,])
      }
      re_union <- re_union[union_genus, union_genus]
      diag(re_union) <- 1
      control_union_data_list <- re_union


      ##sample combination
      case.tax = setdiff(case_dir$X,control_dir$X)
      if (length(case.tax) != 0)
      {
        a = matrix(nrow = length(case.tax),ncol = ncol(control_dir))
        a[,1] = case.tax
        a[,2:ncol(a)] = 0
        a = data.frame(a)
        colnames(a) = colnames(control_dir)
        control_dir.new = rbind(control_dir,a)
      }else
      {
        control_dir.new = control_dir
      }

      control.tax = setdiff(control_dir$X,case_dir$X)
      if (length(control.tax) != 0)
      {
        a = matrix(nrow = length(control.tax),ncol = ncol(case_dir))
        a[,1] = control.tax
        a[,2:ncol(a)] = 0
        a = data.frame(a)
        colnames(a) = colnames(case_dir)
        case_dir.new = rbind(case_dir,a)
      }else
      {
        case_dir.new = case_dir
      }

      if (colnames(case_dir.new)[1] != "X")
      {
        case_dir.new = cbind(rownames(case_dir.new),case_dir.new)
        colnames(case_dir.new)[1] = "X"
      }else
      {
        case_dir.new = case_dir.new
      }

      if (colnames(control_dir.new)[1] != "X")
      {
        control_dir.new = cbind(rownames(control_dir.new),control_dir.new)
        colnames(control_dir.new)[1] = "X"
      }else
      {
        control_dir.new = case_dir.new
      }

      integratedNet = list(case_union_data_list, control_union_data_list, case_dir.new, control_dir.new)

      return(integratedNet)
    } else if (file_test('-f', net_case_dir) &&
               file_test('-f', net_control_dir))
      #####net_case_dir is a file#####
    {
      ###read network
      case_union = read.table(
        net_case_dir,
        header = T,
        sep = '\t',
        row.names = 1
      )
      control_union = read.table(
        net_control_dir,
        header = T,
        sep = '\t',
        row.names = 1
      )

      ##read sample
      case_dir.sample = read.table(
        case_dir,
        header = T,
        sep = '\t'
      )
      control_dir.sample = read.table(
        control_dir,
        header = T,
        sep = '\t'
      )

      ##union
      union_genus1 <- rownames(case_union)
      union_genus2 <- rownames(control_union)

      union_genus <- union(union_genus1, union_genus2)

      union_matrix <-
        matrix(nrow = length(union_genus),
               ncol = length(union_genus),
               0)

      rownames(union_matrix) <- union_genus
      colnames(union_matrix) <- union_genus
      diag(union_matrix) <- 1

      ##combine
      #case
      fi <- case_union
      re_union <- union_matrix
      diff_genus <- setdiff(union_genus, rownames(fi))
      rownames(re_union) <- c(rownames(fi), diff_genus)
      colnames(re_union) <- c(rownames(fi), diff_genus)
      for (j in 1:nrow(fi)) {
        re_union[j, 1:nrow(fi)] <- as.numeric(fi[j,])
      }
      re_union <- re_union[union_genus, union_genus]
      diag(re_union) <- 1
      case_union_data_list <- re_union

      #control
      fi <- control_union
      re_union <- union_matrix
      diff_genus <- setdiff(union_genus, rownames(fi))
      rownames(re_union) <- c(rownames(fi), diff_genus)
      colnames(re_union) <- c(rownames(fi), diff_genus)
      for (j in 1:nrow(fi)) {
        re_union[j, 1:nrow(fi)] <- as.numeric(fi[j,])
      }
      re_union <- re_union[union_genus, union_genus]
      diag(re_union) <- 1
      control_union_data_list <- re_union

      ##sample combination
      case.tax = setdiff(case_dir.sample$X,control_dir.sample$X)
      if (length(case.tax) != 0)
      {
        a = matrix(nrow = length(case.tax),ncol = ncol(control_dir.sample))
        a[,1] = case.tax
        a[,2:ncol(a)] = 0
        a = data.frame(a)
        colnames(a) = colnames(control_dir.sample)
        control_dir.new = rbind(control_dir.sample,a)
      }else
      {
        control_dir.new = control_dir.sample
      }

      control.tax = setdiff(control_dir.sample$X,case_dir.sample$X)
      if (length(control.tax) != 0)
      {
        a = matrix(nrow = length(control.tax),ncol = ncol(case_dir.sample))
        a[,1] = control.tax
        a[,2:ncol(a)] = 0
        a = data.frame(a)
        colnames(a) = colnames(case_dir.sample)
        case_dir.new = rbind(case_dir.sample,a)
      }else
      {
        case_dir.new = case_dir.sample
      }

      if (colnames(case_dir.new)[1] != "X")
      {
        case_dir.new = cbind(rownames(case_dir.new),case_dir.new)
        colnames(case_dir.new)[1] = "X"
      }else
      {
        case_dir.new = case_dir.new
      }

      if (colnames(control_dir.new)[1] != "X")
      {
        control_dir.new = cbind(rownames(control_dir.new),control_dir.new)
        colnames(control_dir.new)[1] = "X"
      }else
      {
        control_dir.new = case_dir.new
      }

      integratedNet = list(case_union_data_list, control_union_data_list, case_dir.new, control_dir.new)

      return(integratedNet)

    } else
      ######else  combination ############
    {
      ######################################1. read data#################################################
      #case data
      setwd(case_dir)
      case_samples <- c()
      case_sample_list <- list()
      dir1 = list.files()
      n1 = length(dir1)
      for (i in 1:n1)
      {
        x.case = read.table(file = dir1[i],
                            header = T,
                            sep = "\t")
        case_samples <-
          c(case_samples, length(colnames(x.case)) - 1)
        case_sample_list[[i]] = x.case
      }

      #control data
      setwd(control_dir)
      control_samples <- c()
      control_sample_list <- list()
      dir2 = list.files()
      n2 = length(dir2)
      for (i in 1:n2)
      {
        x.control = read.table(file = dir2[i],
                               header = T,
                               sep = "\t")
        control_samples <-
          c(control_samples, length(colnames(x.control)) - 1)
        control_sample_list[[i]] = x.control
      }

      #case network
      case_data_list <- list()
      setwd(net_case_dir)
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
      }

      #control network
      control_data_list <- list()
      setwd(net_control_dir)
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
      }


      ##################################union###########################
      #union genus
      union_genus1 <- c()
      for (k in 1:length(case_data_list))
      {
        union_genus1 <- union(union_genus1, rownames(case_data_list[[k]]))
      }

      union_genus2 <- c()
      for (k in 1:length(control_data_list))
      {
        union_genus2 <-
          union(union_genus2, rownames(control_data_list[[k]]))
      }

      union_genus <- union(union_genus1, union_genus2)

      union_matrix <-
        matrix(nrow = length(union_genus),
               ncol = length(union_genus),
               0)

      rownames(union_matrix) <- union_genus
      colnames(union_matrix) <- union_genus
      diag(union_matrix) <- 1

      ##sample combination
      ##case sample
      case.all.sample = data.frame()
      for (k in 1:length(case_sample_list))
      {
        case1 = melt(case_sample_list[[k]],id.vars = "X")
        case.all.sample = rbind(case.all.sample, case1)
      }
      case.com.sample = dcast(case.all.sample,X~variable, mean, fill = 0)

      ##control sample
      control.all.sample = data.frame()
      for (k in 1:length(control_sample_list))
      {
        control1 = melt(control_sample_list[[k]],id.vars = "X")
        control.all.sample = rbind(control.all.sample, control1)
      }
      control.com.sample = dcast(control.all.sample,X~variable, mean, fill = 0)

      case.tax = setdiff(case.com.sample$X,control.com.sample$X)
      if (length(case.tax) != 0)
      {
        a = matrix(nrow = length(case.tax),ncol = ncol(control.com.sample))
        a[,1] = case.tax
        a[,2:ncol(a)] = 0
        a = data.frame(a)
        colnames(a) = colnames(control.com.sample)
        control.com.sample.new = rbind(control.com.sample,a)
      }else
      {
        control.com.sample.new = control.com.sample
      }

      control.tax = setdiff(control.com.sample$X,case.com.sample$X)
      if (length(control.tax) != 0)
      {
        a = matrix(nrow = length(control.tax),ncol = ncol(case.com.sample))
        a[,1] = control.tax
        a[,2:ncol(a)] = 0
        a = data.frame(a)
        colnames(a) = colnames(case.com.sample)
        case.com.sample.new = rbind(case.com.sample,a)
      }else
      {
        case.com.sample.new = case.com.sample
      }


      ######################################2. network construction#################################################
      print ("constructing networks...")
      case_union_data_list <- list()
      for (i in n1) {
        fi <- case_data_list[[i]]
        re_union <- union_matrix
        diff_genus <- setdiff(union_genus, rownames(fi))
        rownames(re_union) <- c(rownames(fi), diff_genus)
        colnames(re_union) <- c(rownames(fi), diff_genus)
        for (j in 1:nrow(fi)) {
          re_union[j, 1:nrow(fi)] <- as.numeric(fi[j,])
        }
        re_union <- re_union[union_genus, union_genus]
        diag(re_union) <- 1
        case_union_data_list[[i]] <- re_union
      }

      control_union_data_list <- list()
      for (i in n2) {
        fi <- control_data_list[[i]]
        re_union <- union_matrix
        diff_genus <- setdiff(union_genus, rownames(fi))
        rownames(re_union) <- c(rownames(fi), diff_genus)
        colnames(re_union) <- c(rownames(fi), diff_genus)
        for (j in 1:nrow(fi)) {
          re_union[j, 1:nrow(fi)] <- as.numeric(fi[j,])
        }
        re_union <- re_union[union_genus, union_genus]
        diag(re_union) <- 1
        control_union_data_list[[i]] <- re_union
      }


      ######################################3. combination#################################################
      print ("integrating networks...")
      pool_union <- union_matrix
      for (i in (1:(nrow(pool_union) - 1))) {
        for (j in ((i + 1):nrow(pool_union))) {
          son <- 0
          mom <- 0
          for (k in n1) {
            v <-
              (1 - case_union_data_list[[k]][i, j] ^ 2) / (case_samples[k] - 1)
            w <- 1 / v
            son <- son + w * case_union_data_list[[k]][i, j]
            mom <- mom + w
          }
          pool_union[i, j] <- son / mom
          pool_union[j, i] <- pool_union[i, j]
        }
      }
      case_union <- pool_union

      pool_union <- union_matrix
      for (i in (1:(nrow(pool_union) - 1))) {
        for (j in ((i + 1):nrow(pool_union))) {
          son <- 0
          mom <- 0
          for (k in n2) {
            v <-
              (1 - control_union_data_list[[k]][i, j] ^ 2) / (control_samples[k] - 1)
            w <- 1 / v
            son <- son + w * control_union_data_list[[k]][i, j]
            mom <- mom + w
          }
          pool_union[i, j] <- son / mom
          pool_union[j, i] <- pool_union[i, j]
        }
      }
      control_union <- pool_union

      if (colnames(case.com.sample.new)[1] != "X")
      {
        case.com.sample.new = cbind(rownames(case.com.sample.new),case.com.sample.new)
        colnames(case.com.sample.new)[1] = "X"
      }else
      {
        case.com.sample.new = case.com.sample.new
      }

      if (colnames(control.com.sample.new)[1] != "X")
      {
        control.com.sample.new = cbind(rownames(control.com.sample.new),control.com.sample.new)
        colnames(control.com.sample.new)[1] = "X"
      }else
      {
        control.com.sample.new = control.com.sample.new
      }

      integratedNet = list(case_union, control_union, case.com.sample.new, control.com.sample.new)

      return(integratedNet)
    }

  }
