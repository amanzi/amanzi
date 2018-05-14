========================================================
How to investigate Continuous Integration testing issues
========================================================
--------------------------------------------------------
1. General Notes
--------------------------------------------------------
- The following website shows the build history for the Amanzi repository
   https://travis-ci.org/amanzi/amanzi/builds

- Login may be required to see it.  If so, make sure your profile on travis-ci.org is tied to the GitHub account you use for Amanzi development.

- There is a time limit of approximately 48 minutes for TravisCI builds
   - This time limit includes both building the code and running all the tests.
   - If a build/test takes too much time, the CI result will show as failed, even if there were no errors.

- Builds that fail with errors will have a red X.

- Builds that fail because they took too long will have a red exclamation point.

- Click on a build to see it's log.

--------------------------------------------------------
2. Particulars of the .travis.yml file
--------------------------------------------------------
- sudo: determines the virtualization environment
   - Setting to false reduces boot time
   - Amanzi's build/test process currently does not need root access

- branches: Specifies which branches of the repository are involved in CI.
   - Amanzi's CI build/tests will only run for commits/pull requests involving the branches specifically listed

- script: the commands to run the tests
   - If any of these commands give a non-zero exit code, the CI fails.
   - Amanzi's script involves the building of a Docker image that both builds the code and runs the tests.
   - This was done to work around the TravisCI time limit
   - It is possible to separate the build stage from the testing stage, however in Amanzi's case this increased the total amount of time spent running the CI such that runs always timed out.

- after_success: Commands run following sucessful build/tests. 
   - This step does not factor into the TravisCI time limit.
   - Settings in TravisCI for the amanzi/amanzi repository include setting the values for:
      - DOCKERHUB
      - DOCKER_PASSWORD
      - DOCKER_USERNAME
      - DEPLOY
   - TravisCI is currently set to push the newly-built Amanzi Docker image to the appropriate DockerHub repository.
   - This ensures that the DockerHub image tagged metsi/amanzi:latest truly is based on the latest commit to the master branch.
