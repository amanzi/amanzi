  // first figure out how many entries there are
  int nblocks = 0;
  for (Teuchos::ParameterList::ConstIterator i = wrm_plist.begin(); i != wrm_plist.end(); i++) {
    // only count sublists
    if (wrm_plist.isSublist(wrm_plist.name(i))) {
	  nblocks++;
	} else {
	  // currently we only support van Genuchten, if a user
	  // specifies something else, throw a meaningful error...
	  Errors::Message m("PermafrostProblem: the Water retention models sublist contains an entry that is not a sublist!");
	  Exceptions::amanzi_throw(m);
	}
  }

  WRM_.resize(nblocks);

  int iblock = 0;
  for (Teuchos::ParameterList::ConstIterator i = wrm_plist.begin(); i != wrm_plist.end(); i++) {
    if (wrm_plist.isSublist(wrm_plist.name(i))) {
      Teuchos::ParameterList &block_wrm_plist = wrm_plist.sublist( wrm_plist.name(i) );

      // which water retention model are we using? currently we only have van Genuchten
      if ( block_wrm_plist.get<string>("Water retention model") == "van Genuchten") {
        // read the mesh block number that this model applies to
        int meshblock = block_wrm_plist.get<int>("Region ID");

        // read values for the van Genuchten model
        double vG_m = block_wrm_plist.get<double>("van Genuchten m");
        double vG_alpha = block_wrm_plist.get<double>("van Genuchten alpha");
        double vG_sr = block_wrm_plist.get<double>("van Genuchten residual saturation");

        WRM_[iblock] = Teuchos::rcp(new vanGenuchtenModel(meshblock,vG_m,vG_alpha,
                                                          vG_sr,p_atm_));
      }
      iblock++;
    }
  }
