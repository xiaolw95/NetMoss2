#' @title Identifying Differential Bacteria Based on a Network Workflow.
#' @description This is a function to identify differential bacteria between case and control data set.
#'
#' @param case_dir string.The directory of diseased data set or a single diseased file.
#' @param control_dir string.The directory of healthy data set or a single healthy file.
#' @param net_case_dir string.The directory of network correlation of diseased data set or a single diseased network file.
#' @param net_control_dir string.The directory of network correlation of healthy data set or a single healthy network file.
#' @param scaled logical.If TURE then the NetMoss score are scaled.
#' @param deepSplit numerical. For method "hybrid", can be either logical or integer in the range 0 to 4. For method "tree", must be logical. In both cases, provides a rough control over sensitivity to cluster splitting. The higher the value (or if TRUE), the more and smaller clusters will be produced. For the "hybrid" method, a finer control can be achieved via maxCoreScatter and minGap below.
#' @param minModuleSize numerical. Minimum module size. Set to 20 by default.
#' @param soft_pow numerical. Soft power for module division. Set to 5 by default.
#'
#' @return NetMoss score and module division results of each bacterium.
#' @export
#'
#' @examples
#' data(testData)
#' nodes_result = NetMoss(case_dir = mydata[[1]],
#'      control_dir = mydata[[2]],
#'      net_case_dir = mydata[[3]],
#'      net_control_dir = mydata[[4]])
#'
NetMoss <-
  function(case_dir,
           control_dir,
           net_case_dir,
           net_control_dir,
           scaled = TRUE,
           deepSplit = 4,
           minModuleSize = 20,
           soft_pow = 5) {
    my.wd = getwd()

    result = list()

    netAll = getNetwork(
      case_dir = case_dir,
      control_dir = control_dir,
      net_case_dir = net_case_dir,
      net_control_dir = net_control_dir
    )

    case_union = netAll[[1]]
    control_union = netAll[[2]]
    case_sample_union = netAll[[3]]
    control_sample_union = netAll[[4]]

    modAll = divModule(case_union = case_union,
                       control_union = control_union,
                       minModuleSize = minModuleSize,
                       deepSplit = deepSplit,
                       soft_pow = soft_pow)
    dissTOM_case = modAll[[1]]
    dissTOM_control = modAll[[2]]
    dynamicMods_case = modAll[[3]]
    dynamicMods_control = modAll[[4]]

    nodes_result = NetzGO(
      control_mat = control_union,
      case_mat = case_union,
      control_dist = dissTOM_control,
      case_dist = dissTOM_case,
      control_mod = as.numeric(dynamicMods_control),
      case_mod = as.numeric(dynamicMods_case),
      scaled = scaled
    )

    setwd(my.wd)

    result[[1]] = nodes_result
    result[[2]] = case_union
    result[[3]] = control_union
    result[[4]] = case_sample_union
    result[[5]] = control_sample_union

    message ("NetMoss analysis completed.")

    return(result)

  }
