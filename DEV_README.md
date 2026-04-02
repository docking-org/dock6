# dock6

Trent Balius (document based on DOCK3.7 by Ryan G. Coleman Teague Sterling)

DOCK6 is a molcular docking program written in C++.  

How To Fork for Development
---------------------------
 0. Create an account on github.com

 1. Go to https://github.com/settings/ssh and paste your public key in the box

 2. Go to https://github.com/docking-org/DOCK

 3. Click "Fork" to create your own working copy

 4. On your development machine and set the the trunk as upstream:

       git clone git@github.com:[YOUR_USERNAME_HERE]/dock6.git dock6

       cd dock6

       git remote add upstream git@github.com:docking-org/dock6.git

 5. After making some changes, save them locally on your machine:

       git commit -m "I changed test/file-abc" test/file-abc

    or

       git commit -a -m "I changed all the files"

 6. After a series of changes are ready to push remotely

       git push

 7. If these changes should be added to the master DOCK branch

    - Go to https://github.com:[YOUR_USERNAME_HERE]/dock6.git
    - Click "Pull Request" in the upper-right corner of the file listing
    - Click "Click to create a pull request for this comparison"
    - Add any additional notes you may want
    - Press "Send pull request"

 8. If you are an administrator in the primary DOCK repository, click "Merge pull request" , add any comments wanted, then click "Confirm merge" and the changes will be added to the DOCK trunk

 9. To update your fork with the trunk

       git pull upstream master
       git push

Directories
-----------

    test: Tools and data for testing dockenv functionality

Usage
-----

**Note:** This is just a stub of the usage documentation. Consult the wiki linked above for more information.

You will need to set the environmental variable DOCKBASE to the root repository path.

For example usage look the test directory.
