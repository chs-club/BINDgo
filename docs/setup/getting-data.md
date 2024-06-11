# Downloading the data

In this tutorial, we will explore how to download specific files from a Kaggle dataset using the Kaggle API.
The Kaggle API provides a convenient way to access and retrieve datasets programmatically, allowing you to automate the process and integrate it into your data pipeline.
We will use the dataset located at [https://www.kaggle.com/competitions/leash-BELKA](https://www.kaggle.com/competitions/leash-BELKA) as an example and focus on downloading only the `train.parquet` and `test.parquet` files.

## Prerequisites

First and foremost, you need to have a Kaggle account.
If you don't have one already, head over to the [Kaggle website](https://www.kaggle.com/) and create an account.

### Kaggle API

To interact with the Kaggle API, we'll be using the `kaggle` package.
You can install this package by running the following command in your terminal or command prompt:

```bash
pip install kaggle
```

This command will download and install the Kaggle API package, along with its dependencies, making it ready for use.
You can find more information in the [documentation](https://github.com/Kaggle/kaggle-api/blob/main/docs/README.md).

### Kaggle API key

Once you have an account, you'll need to obtain your API credentials, which consist of a username and a key.
These credentials are essential for authenticating your requests to the Kaggle API.
To do this, click on your Kaggle profile photo on the upper right and select `Settings`.
Under `Account`, there is an `API` section that will allow you to `Create New Token`.
This will download a file called `kaggle.json`.

```json
{"username":"aalexmmaldonado","key":"1c8b25b85741554a1edd12a7c0c608de"}
```

!!! danger

    It's important to keep this file secure and not share it with anyone, as it grants access to your Kaggle account.
    The one above is a randomly generated token just for demonstrative purposes.

Next, you need to place the `kaggle.json` file in a specific location on your machine.
For Windows users, the file should be placed at `C:\Users\<username>\.kaggle\kaggle.json`.
For Linux and macOS users, the file should be located at `~/.kaggle/kaggle.json`.
Make sure to replace `<username>` with your actual username on your machine.

!!! note

    You often have to create the directory `.kaggle` if it does not already exist.

To ensure that the `kaggle.json` file has the appropriate permissions, open your terminal or command prompt and run the following command:

```bash
chmod 600 ~/.kaggle/kaggle.json
```

This command sets the permissions of the `kaggle.json` file to be readable and writable only by the owner, providing an extra layer of security.

## Downloading

With the API credentials set up, we're now ready to download the specific files we need.
Open your terminal or command prompt and navigate to the directory where you want to download the files.
This can be done using the `cd` command followed by the desired directory path.

To download the `test.parquet.zip` file (18.2 MB), run the following command:

```bash
kaggle competitions download -c leash-BELKA -f test.parquet && unzip test.parquet.zip
```

Similarly, to download the `train.parquet` file (1.62 GB), run:

```bash
kaggle competitions download -c leash-BELKA -f train.parquet
```

In these commands, the `-c` flag specifies the competition or dataset name, and the `-f` flag specifies the file name you want to download. The Kaggle API will authenticate your request using the credentials from the `kaggle.json` file and initiate the download process.

Wait for the download to complete. The files will be saved in the current directory. Once the download is finished, you can verify that the `test.parquet` and `train.parquet` files have been successfully downloaded by checking the contents of the current directory.

With the specific files now available locally, you can proceed to use them in your data analysis or machine learning projects. Load the files into your preferred data processing or analysis tools, such as pandas, and start exploring and extracting insights from the data.

## Shortcut

If you are using this repository, you can run

```bash
make get-data
```

and it will download and extract the data for you in the `data/` directory.
