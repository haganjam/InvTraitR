#
# The tests run in isolation, but we need our database files
# for the tests to sucessfully execute. Hence, we need to fetch
# them for the tests. Downloading them every time a test execute
# would be super wasteful - so let's grab them.
#
download_database(appdata_subdir)
