#' Convenient computation of trigonometric predictor
#'
#' @param coef Matrix of beta coefficients
#' @param time timepoint to evaluate
#' @param degree number of trigonometric functions to use
#' @param L length of one period
#'
#' @return
#' @export
#'
#' @examples
pv = function(coef, time, degree, L = 24){
  if(degree == 1){
    return(coef[,1] + coef[,2]*sin(2*pi*time/L) + coef[,3]*cos(2*pi*time/L))
  } else if(degree == 2){
    return(coef[,1] + coef[,2]*sin(2*pi*time/L) + coef[,3]*sin(2*pi*time*2/L) +
             coef[,4]*cos(2*pi*time/L) + coef[,5]*cos(2*pi*time*2/L))
  } else if(degree == 3){
    return(coef[,1] + coef[,2]*sin(2*pi*time/L) + coef[,3]*sin(2*pi*time*2/L) + coef[,4]*sin(2*pi*time*3/L) +
             coef[,5]*cos(2*pi*time/L) + coef[,6]*cos(2*pi*time*2/L) + coef[,7]*cos(2*pi*time*3/L))
  } else if(degree == 4){
    return(coef[,1] + coef[,2]*sin(2*pi*time/L) + coef[,3]*sin(2*pi*time*2/L) + coef[,4]*sin(2*pi*time*3/L) + coef[,5]*sin(2*pi*time*4/L) +
             coef[,6]*cos(2*pi*time/L) + coef[,7]*cos(2*pi*time*2/L) + coef[,8]*cos(2*pi*time*3/L) + coef[,9]*cos(2*pi*time*4/L))
  } else if(degree == 5){
    return(coef[,1] + coef[,2]*sin(2*pi*time/L) + coef[,3]*sin(2*pi*time*2/L) + coef[,4]*sin(2*pi*time*3/L) + coef[,5]*sin(2*pi*time*4/L) + coef[,6]*sin(2*pi*time*5/L) +
             coef[,7]*cos(2*pi*time/L) + coef[,8]*cos(2*pi*time*2/L) + coef[,9]*cos(2*pi*time*3/L) + coef[,10]*cos(2*pi*time*4/L) + coef[,11]*cos(2*pi*time*5/L))
  }
}
