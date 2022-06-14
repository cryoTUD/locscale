import unittest

def download_test_data_from_url(download_folder):
    import wget
    print("Downloading test data \n")
    #url_test_data = "https://surfdrive.surf.nl/files/index.php/s/xJKxGXR0LWGBDWM/download"
    url_test_data = "https://surfdrive.surf.nl/files/index.php/s/lk9CdNO5gszFll1/download"
    
    wget.download(url_test_data, download_folder)

def extract_tar_files_in_folder(tar_folder, use_same_folder=True):
    import tarfile
    import os
    if use_same_folder == 0:
        target_folder = tar_folder
    else:
        target_folder = os.path.dirname(tar_folder)

    for file in os.listdir(tar_folder):
        if file.endswith(".tar.gz"):
            print("Extracting: {}".format(file))
            tar = tarfile.open(os.path.join(tar_folder,file))
            tar.extractall(target_folder)
            tar.close()

def download_and_test_everything():

    import tempfile
    import unittest
    import os
    import locscale
    ## Create folder to download tests_data
    test_data_path = os.path.join(os.path.dirname(locscale.__file__), "tests","test_data")
    test_path = os.path.join(os.path.dirname(locscale.__file__), "tests")
    if not os.path.exists(test_data_path):
        os.makedirs(test_data_path, exist_ok=True)
    
        ## Download test data
        download_test_data_from_url(test_data_path)
        ## Extract tar files
        extract_tar_files_in_folder(test_data_path)

    ## Create test suite
    test_loader = unittest.TestLoader()
    test_suite = test_loader.discover(test_path, pattern="test_*.py")
    ## Run test suite
    test_runner = unittest.TextTestRunner(verbosity=2)
    test_result = test_runner.run(test_suite)
    ## Print test result
    print("Test result: {}".format(test_result))
    if test_result.wasSuccessful():
        print("All tests passed")
        
    else:
        print("Some tests failed")
        print(test_result.printErrors())
    
    print("="*80)
    print("Done")
    print("="*80)




