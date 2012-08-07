==============================================
Amanzi Verification and Validation
==============================================

Borrow some overview stuff from the Basic Phase Test Plan document

Testing Framework  
-------------------

The CTest/CDash tools provide a basis for a creating a flexible testing framework that supports a hierarchical test suite with
automated execution and reporting capabilities.  The description of the testing framework, test specification and usage
is described ... (not sure we need this here)


Basic Phase Testing 
---------------------

We are collecting tests from the literature to test flow, transport, and chemistry.  The organisation of these tests is still 
under consideration, but for now they are nominally organised by the dominant process of interest:

   * Flow (steady-state and transient single-phase saturated and unsaturated models)
   * Transport
   * Chemistry

These tests are run as part of our nightly (soon I hope) testing on multiple platforms. But they also provide a good suite
of simple examples for you to run and experiment with.  

