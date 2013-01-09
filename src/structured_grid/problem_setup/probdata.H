      REAL_T  denerr, adverr
      REAL_T  rhoval(1:100)
      REAL_T  muval(1:100)
      REAL_T  gravity
      integer phaseidx(1:100)
      integer ncomps, nphases, nreact, ntracers, naux
      integer probtype, model

      common /probdt1/ probtype, model
      common /probdt2/ denerr, adverr
      common /probdt3/ ncomps, nphases, nreact, ntracers, naux
      common /probdt4/ rhoval, muval, phaseidx
      common /probdt5/ gravity

      namelist /fortin/ probtype, denerr, adverr
      save /probdt1/, /probdt2/, /probdt3/, /probdt4/, /probdt5/

      REAL_T sThresh
      integer use_shifted_kr
      common /eoscomr/ sThresh
      common /eoscomi/ use_shifted_kr
      save /eoscomr/, /eoscomi/
