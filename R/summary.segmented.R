`summary.segmented` <-
function(object, short=FALSE, ...){
#revisione 13/05/03;7/10/03;28/11;24/02/04
    if(is.null(object$psi)) object<-object[[length(object)]]
#i seguenti per calcolare aa,bb,cc funzionano per lm e glm, da verificare con arima....
#    nome<-rownames(object$psi)
#    nome<-as.character(parse("",text=nome))
#    aa<-grep("U",names(coef(object)[!is.na(coef(object))]))
#    bb<-unlist(sapply(nome,function(x){grep(x,names(coef(object)[!is.na(coef(object))]))},simplify=FALSE,USE.NAMES=FALSE))
#    cc<-intersect(aa,bb) #indices of diff-slope parameters
#    iV<- -grep("psi.",names(coef(object)[!is.na(coef(object))]))#indices of all but the Vs
    nomiU<-object$nameUV[[1]]
    idU<-match(nomiU,names(coef(object)[!is.na(coef(object))]))
    nomiV<-object$nameUV[[2]]
    idU<-match(nomiU,names(coef(object)[!is.na(coef(object))]))
    idV<-match(nomiV,names(coef(object)[!is.na(coef(object))]))
    beta.c<- coef(object)[nomiU]
    if("lm"%in%class(object) && !"glm"%in%class(object)){
        summ <- c(summary.lm(object, ...), object["psi"])
        summ$Ttable<-summ$coefficients[-idV,] #ho tolto -4 per inserire anche i p-value
        summ$Ttable[idU,4]<-NA
        coeff<-summ$coefficients[,1]
        v<-summ$coefficients[,2]
        summ$gap<-cbind(coeff[idV]*beta.c,v[idV]*beta.c,coeff[idV]/v[idV])
        #dimnames(summ$gap)<-list(rep("",nrow(object$psi)),c("Est.","SE","t value"))
        colnames(summ$gap)<-c("Est.","SE","t value")
        rownames(summ$gap)<-nomiU
        summ[c("it","epsilon")]<-object[c("it","epsilon")]
        summ$short<-short
        class(summ) <- c("summary.segmented", "summary.lm")
        return(summ)}
    if("glm"%in%class(object)){
        summ <- c(summary.glm(object, ...), object["psi"])
        summ$Ttable<-summ$coefficients[-idV,]
        summ$Ttable[idU,4]<-NA
        coeff<-summ$coefficients[,1]
        v<-summ$coefficients[,2]
        summ$gap<-cbind(coeff[idV]*beta.c,v[idV]*beta.c,coeff[idV]/v[idV])
        #dimnames(summ$gap)<-list(rep("",nrow(object$psi)),c("Est.","SE","t value"))
        colnames(summ$gap)<-c("Est.","SE","t value")
        rownames(summ$gap)<-nomiU
        summ[c("it","epsilon")]<-object[c("it","epsilon")]
        summ$short<-short
        class(summ) <- c("summary.segmented", "summary.glm")
        return(summ)}
    if("Arima"%in%class(object)){
        #da controllare
        coeff<-object$coef
        v<-sqrt(diag(object$var.coef))
        Ttable<-cbind(coeff[-idV],v[-idV],coeff[-idV]/v[-idV])
        object$gap<-cbind(coeff[idV]*beta.c,v[idV]*beta.c,coeff[idV]/v[idV])
        #dimnames(object$gap)<-list(rep("",nrow(object$psi)),c("Est.","SE","t value"))
        colnames(summ$gap)<-c("Est.","SE","t value")
        rownames(summ$gap)<-nomiU
        colnames(Ttable)<-c("Estimate","Std. Error","t value")
        object$Ttable<-Ttable
        object$short<-short
        summ<-object 
        class(summ) <- "summary.segmented"
        return(summ)}
}

