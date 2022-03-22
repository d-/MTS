RcppExport SEXP GetVarmaResiduals(SEXP _timeSeries,
                                  SEXP _mask,
                                  SEXP _paramFixed,
                                  SEXP _p,
                                  SEXP _q,
                                  SEXP _isMeanIncluded
                                  );

RcppExport SEXP GetSVarmaResiduals(SEXP _timeSeries,
                                   SEXP _mask,
                                   SEXP _paramFixed,
                                   SEXP _orders,
                                   SEXP _arlags,
                                   SEXP _malags,
                                   SEXP _sresi,
                                   SEXP _swi,
                                   SEXP _isMeanIncluded
                                   );

RcppExport SEXP GetVMAObs(SEXP _timeSeries,
                          SEXP _mask,
                          SEXP _paramFixed,
                          SEXP _q,
                          SEXP _isMeanIncluded
                          );

RcppExport SEXP GetVMATH(SEXP _timeSeries,
                         SEXP _mask,
                         SEXP _paramFixed,
                         SEXP _q,
                         SEXP _isMeanIncluded
                         );
