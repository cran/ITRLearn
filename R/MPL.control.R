MPL.control <-
function(pi.est = NULL, h.est = NULL, boot.sample = 600)
{
    rval <- list(pi.est=pi.est, h.est=h.est, boot.sample=boot.sample)
    rval
}