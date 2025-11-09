# GitLab # {#dvGitlab}

GitLab provides an entire ecosystem for the collaborative development of software using the version control protocol __git__.
Here, we focus on the relevant features for the work flow used for the development of @maia.

## Issues
If you want to report some unwanted behavior of the code __or__ propose ideas for futures development ideas, use the __issue__ function.
Please provide a concise description of the problem or your idea. You should use our predefined labels to categorize your issue. There are labels referring to the different parts of the code as well as labels expressing the point of the issue.
Feel free to use the issue as a central place to organize collaborative efforts by assigning other members and exchange your ideas in the comment section.

After you have created an issue you have the option to directly create a branch for this specific issue in the web-interface. This is highly recommended, since the branch name is set automatically including the number of the issue.
This makes life much easier down the line, because GitLab can assign your issue to the branch (and later to the merge request). 

## Merge Requests

If you want to propose your changes to the code to be included into the master, you have to open a merge request (MR).
You should provide a comprehensive description for you merge request, stating which parts of the code will be affected by your changes and in what way. Again, use the predefined labels to categorize your proposed changes.
Other members can now comment on your changes. Please answer to the comments accordingly and make changes to your code if needed.
As soon as your code is approved, it will be merged into the master branch of the code. Your branch will be deleted and assigned issues will be marked as resolved. 

## Continuous integration

To ensure the code is running as expected using all relevant compilers and configurations, we make use of the continuous integration features of GitLab.
On different occasions we run a testing pipeline to check the validity of the code. This pipeline consists of different staged which in turn a composed of multiple jobs.
While being triggered by the GitLab instance, the job itself is executed by connected GitLab runner instances at the Institute of Aerodynamics.

The most important stages/jobs fall into one of the following categories:

- __Build jobs:__ These jobs do nothing but check if the code can be compiled using different compilers with different build types. Additionally we remotely test the compilation on selected HPC machines.
- __Test jobs:__ Using our [testing environment](@ref dvTesting), the validity of the code is checked for different compilers/build types.
- __Deploy jobs:__ Packages like the ParaView plugins depend on @maia to be build. That is why the code is deployed there to run the respective CI operations. Additionally, this documentation is deployed.

### Merge request pipeline
When pushing to a branch for which an open MR exists, a small pipeline with a reduced number of building and testing jobs is started. If your branch is still under development, you don't want this pipeline to run every time you push some changes, because the computational resource needed for the execution are not for free and could be used in a more meaningful way. To prevent the pipeline from running you can mark your MR as work-in-progress by prepending "Draft:" or "WIP:" to the title. If you are finished with your changes, just remove this tag again.

@note Unfortunately, the pipeline does not start after you mark your MR as ready. You have to push a new (dummy) commit to your branch to trigger the pipeline.

A successful MR pipeline is a necessary criterion for your changes to be merged to the master.

### Master merge pipeline

After a MR has been tested, reviewed and approved, the branch is merged to the master. This triggers an extensive pipeline testing many configurations. This pipeline should not fail if possible. If, for some unforeseeable reason, your pipeline fails, you are responsible for fixing this __as soon as possible__. Do this by pushing the fixes to another branch, ideally named something like **fix_old_branch_name**, and open a new MR.
