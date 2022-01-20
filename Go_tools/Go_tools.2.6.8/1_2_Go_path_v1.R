#' A Go_huamnn2ps
#'
#'
#' @param huamnn2ps
#' @keywords huamnn2ps
#' @export
#' @examples
#' Go_huamnn2ps


Go_path <- function(project, pdf, table, path){
  # main dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  

  # main pdf
  if(is.null(pdf)){
    print("No pdf dir.")
  } else if (pdf == "yes" | pdf == "Yes"|pdf == "YES"){
    out_pdf <- file.path(sprintf("%s/pdf",out)) 
    if(!file_test("-d", out_pdf)) dir.create(out_pdf)
    
    print("pdf is in your working dir. Use /dir$pdf/ for save.")
  }

  
  
  # main table
  if(is.null(table)){
    print("No table dir.")
  } else if (table == "yes" | table == "Yes"|table == "YES"){
    out_tab <- file.path(sprintf("%s/table",out)) 
    if(!file_test("-d", out_tab)) dir.create(out_tab)

    print("table is in your working dir.Use /dir$tab/ for save.")
  }
  
  if(is.null(path)){
    print("No another dir.")
  } else if(!is.null(path)){
    out_path <- file.path(sprintf("%s/%s",out,path)) 
    if(!file_test("-d", out_path)) dir.create(out_path)
    print("path is in your working dir. Use /dir$path/ for save.")
  }


  # 한개 이상 return 하기
  
  
  functionReturningTwoValues <- function() {
    dirs <- list()
    if(is.null(pdf)){
    } else if (pdf == "yes" | pdf == "Yes"|pdf == "YES"){
      dirs$pdf <- out_pdf
    }
    if(is.null(table)){
      next
    } else if (table == "yes" | table == "Yes"|table == "YES"){
      dirs$tab <- out_tab
    }
    if(is.null(path)){
    } else if(!is.null(path)){
      dirs$path <- out_path
    }

    return(dirs) 
  }

  functionReturningTwoValues ()
  
}
