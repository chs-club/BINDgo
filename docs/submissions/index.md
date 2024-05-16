# Submissions

TODO:

## Starting your submission

This guide will walk you through the steps to create your own submission using Git branches and Jupyter notebooks.
By following these instructions, you'll be able to develop and submit your machine learning model independently.

Prerequisites:

-   Git installed on your local machine.
-   Access to the challenge repository.
-   Jupyter Notebook environment set up.

Notes:

-   Keep your branch separate from the main branch and other students' branches. This allows everyone to work independently on their own submissions.
-   Regularly commit and push your changes to ensure your work is backed up and accessible.
-   If you encounter any conflicts or issues, reach out to the challenge organizers or seek assistance from your peers.

### Clone the Challenge Repository

1.  Open your terminal or Git client.
2.  Navigate to the directory where you want to clone the challenge repository.
3.  Run the following command to clone the repository:

   ```bash
   git clone https://github.com/oasci/leash-bio-kaggle
   ```

### Create a New Branch

1.  Navigate to the cloned repository directory:

   ```bash
   cd leash-bio-kaggle
   ```

   Replace `<repository-directory>` with the name of the cloned repository directory.
2. Create a new branch for your submission using the following command:

   ```bash
   git checkout -b your_name
   ```

   Replace `your-name` with your actual name or a unique identifier.
   This branch will be used to track your own submission separately from the main branch and other participants' branches.

### Copy the Example Submission Directory

1.  In the repository, locate the directory named `jane_doe`.
    This directory serves as an example submission structure.
2.  Create a copy of this directory and rename it with your own name. Run the following command:

    ```bash
    cp -r jane_doe your_name
    ```

    Replace `your_name` with the same name you used for your branch in Step 2.

### Develop Your Machine Learning Model

1.  Navigate to your newly created submission directory:

    ```bash
    cd your-name
    ```
2.  Create a new Jupyter notebook for your machine learning model:

    ```bash
    jupyter notebook
    ```
    This command will open Jupyter Notebook in your web browser.
3.  In the Jupyter Notebook interface, click on "New" and select "Python 3" to create a new notebook.
4.  Rename the notebook to a descriptive name for your model.
5.  Develop and train your machine learning model within the notebook. You can create additional notebooks as needed.
6.  Save your notebooks regularly to ensure your work is preserved.

### Commit and Push Your Changes

1.  After making changes to your notebooks or adding new files, commit your changes using the following commands:

    ```bash
    git add .
    git commit -m "Add your commit message here"
    ```

    Replace `"Add your commit message here"` with a descriptive message summarizing your changes.
2.  Push your changes to the remote repository:

    ```bash
    git push origin your-name
    ```

    Replace `your-name` with the name of your branch.

### Iterate and Update

1.  Continue working on your notebooks, making improvements and iterations to your machine learning model.
2.  Regularly commit and push your changes to keep your branch up to date.
3.  If you need to switch between your branch and the main branch, use the following commands:

    ```bash
    git checkout main  # Switch to the main branch
    git checkout your-name  # Switch back to your branch
    ```
