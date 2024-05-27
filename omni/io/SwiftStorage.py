"""Swift class for remote storage."""
from swiftclient.service import SwiftService, SwiftError, Connection
from swiftclient.multithreading import MultiThreadingManager
import datetime
import dateutil.parser
import requests
from bs4 import BeautifulSoup
import re
import logging
from omni.io.RemoteStorage import RemoteStorage
logging.basicConfig(level=logging.ERROR)
logging.getLogger("requests").setLevel(logging.DEBUG)
logging.getLogger("swiftclient").setLevel(logging.DEBUG)
logger = logging.getLogger(__name__)


def get_swift_auth(user, key, authurl, tenant_name, auth_version):
    """
    Get the Swift authentication url and token from user password.

    Args:
        user (str): The username.
        key (str): The user password.
        authurl (str): The swift URL for authentication.
        tenant_name (str): The name of the tenant or project.
        auth_version (str): The version of the authentication service.

    Returns:
        dict: A dictionary containing the Swift authentication url and token.
    """
    # temporary connection to get auth token
    tmp_con = Connection(authurl=authurl, user=user, key=key, tenant_name=tenant_name, auth_version=auth_version)
    auth = tmp_con.get_auth()
    swift_auth_options = {'preauthurl': auth[0], 'preauthtoken': auth[1], 'retries': 3, 'starting_backoff': 1, 'max_backoff':8}
    return swift_auth_options


def get_meta_mtime(preauthurl, containername, objectname):
    """
    Retrieves the metadata modification time and access time of a file from a Swift storage.

    Args:
        preauthurl (str): The pre-authenticated URL of the Swift storage.
        containername (str): The name of the container where the file is stored.
        objectname (str): The name of the file.

    Returns:
        tuple: A tuple containing the file's metadata modification time and access time.

    Raises:
        requests.HTTPError: If the HTTP request to retrieve the file fails.
    """
    urlfile = f"{preauthurl}/{containername}/{objectname}"
    response = requests.get(urlfile)
    if response.ok:
        response_headers = response.headers
        if 'X-Object-Meta-Mtime' in response_headers.keys():
            file_time = datetime.datetime.fromtimestamp(float(response_headers['X-Object-Meta-Mtime']), datetime.timezone.utc)
        elif 'X-Object-Meta-Last-Modified' in response_headers.keys():
            file_time = dateutil.parser.parse(response_headers['X-Object-Meta-Last-Modified'])
        else:
            file_time = dateutil.parser.parse(response_headers['Last-Modified'])
        access_time = dateutil.parser.parse(response_headers['Date'])
        return file_time, access_time
    else:
        response.raise_for_status()



class SwiftStorage(RemoteStorage):
    def __init__(self, auth_options, benchmark):
        super().__init__(auth_options, benchmark)
        self.containers = list()
        if 'preauthtoken' in self.auth_options.keys():
            self._test_connect()
            self._get_benchmarks()

            if not benchmark in self.benchmarks:
                logger.warning(f"Benchmark {benchmark} does not exist, creating new benchmark.")
                self._create_benchmark(benchmark)
                self._get_containers()
                self._get_benchmarks()

            self._get_versions()
        else:
            self._get_versions(readonly=True)
            # read-only mode
            pass
    

    def connect(self):
        """
        Connects to the Swift storage.

        Returns:
            SwiftService: The SwiftService object representing the connection to the Swift storage.
        """
        if 'preauthurl' in self.auth_options.keys() and 'preauthtoken' in self.auth_options.keys():
            # from: https://opendev.org/openstack/python-swiftclient/src/branch/master/swiftclient/service.py
            swift = SwiftService()
            def create_connection():
                return Connection(**self.auth_options)
            swift.thread_manager = MultiThreadingManager(
                        create_connection,
                        segment_threads=swift._options['segment_threads'],
                        object_dd_threads=swift._options['object_dd_threads'],
                        object_uu_threads=swift._options['object_uu_threads'],
                        container_threads=swift._options['container_threads']
                    )
            return swift
        else:
            return SwiftService(options=self.auth_options)
    
    def _test_connect(self):
        with self.connect() as swift:
            stat = swift.stat()
            if not stat['success']:
                raise SwiftError(f"Connection failed: {stat['error']}")



    def _get_containers(self):
        """
        Retrieves the list of containers from the Swift storage.

        Returns:
        - A list of container names.
        """
        with self.connect() as swift:
            try:
                containers = list()
                for page in swift.list():
                    for element in page['listing']:
                        containers.append(element['name'])
                self.containers = containers
                return self.containers

            except SwiftError as e:
                logger.error(e.value)

   
    def _get_benchmarks(self, update=True):
        if len(self.containers) == 0 or update:
            self._get_containers()
        benchmarks = list()
        for con in self.containers:
            # remove major and minor version from benchmark name
            bm = ''.join(con.split(".")[:-2])
            if bm not in benchmarks and bm != "":
                benchmarks.append(bm)
        self.benchmarks = benchmarks
        return self.benchmarks 
 
    def _create_benchmark(self, benchmark, update=True):
        if len(self.benchmarks) == 0 or update:
            self._get_benchmarks()
        if benchmark in self.benchmarks:
            raise ValueError("Benchmark already exists")
        # create new version
        with self.connect() as swift:
            post = swift.post(container=f"{benchmark}.test.1", options={'read_acl': ".r:*,.rlistings"})
            if not post['success']:
                raise SwiftError(f"Benchmark creation of {benchmark}.test.1 failed: {post['error']}")
            post = swift.post(container=f"{benchmark}.overview", options={'read_acl': ".r:*,.rlistings"})
            if not post['success']:
                raise SwiftError(f"Benchmark creation of {benchmark}.overview failed: {post['error']}")
            post = swift.post(container=f"{benchmark}.0.1", options={'read_acl': ".r:*,.rlistings"})
            if not post['success']:
                raise SwiftError(f"Benchmark creation of {benchmark}.0.1 failed: {post['error']}")

        self._update_overview(cleanup=True)
            # except SwiftError as e:
            #     logger.error(e.value)
        
    def _get_versions(self, update=True, readonly=False):
        if not 'preauthtoken' in self.auth_options.keys() or readonly:
            url = f"{self.auth_options['preauthurl']}/{self.benchmark}.overview"
            params={'format': 'xml'}
            response = requests.get(url, params=params)
            if response.ok:
                response_text = response.text
            else:
                response.raise_for_status()
            soup = BeautifulSoup(response_text, 'xml')
            self.versions = [obj.find('name').text for obj in soup.find_all('object')]
            return self.versions

            # with self.connect() as swift:
            #     versions = list()
            #     for page in swift.list(container=f"{self.benchmark}.overview"):
            #         for element in page['listing']:
            #             versions.append(element['name'])
            # self.versions = versions
            # return self.versions
        else:
            if update:
                self._get_containers()
            versions = list()
            for con in self.containers:
                if re.match(f"{self.benchmark}.(\\d+).(\\d+)", con):
                    versions.append('.'.join(con.split(".")[-2:]))
            self.versions = versions
            return self.versions
    
    def _update_overview(self, cleanup=True):
        # check versions in overview
        versions_in_overview = []
        with self.connect() as swift:
            for page in swift.list(container=f"{self.benchmark}.overview", options={'long' :  True}):
                for element in page['listing']:
                    versions_in_overview.append(element['name'])
        # available benchmark versions
        self._get_versions()
        # missing versions
        versions_not_in_overview = [v for v in self.versions if v not in versions_in_overview ]
        # create missing versions
        with self.connect() as swift:
            for v in versions_not_in_overview:
                post = swift.post(container=f"{self.benchmark}.overview/{v}")
                if not post['success']:
                    raise SwiftError(f"Object creation of {self.benchmark}.overview/{v} failed: {post['error']}")

        # remove unavailable versions 
        if cleanup:
            versions_in_overview = []
            with self.connect() as swift:
                for page in swift.list(container=f"{self.benchmark}.overview", options={'long' :  True}):
                    for element in page['listing']:
                        versions_in_overview.append(element['name'])

            versions_in_overview_but_unavailable = [v for v in versions_in_overview if v not in self.versions ]
            with self.connect() as swift:
                deletes = swift.delete(container=f"{self.benchmark}.overview",objects=versions_in_overview_but_unavailable)
                for delete in deletes:
                    if not delete['success']:
                        raise SwiftError(f"Object deletion of {self.benchmark}.overview/{v} failed: {delete['error']}")


   
        
    def _create_new_version(self):
        if self.major_version_new is None or self.minor_version_new is None:
            raise ValueError("No version provided")
        
        # update
        self._get_versions()
        # check if version exists
        if not f"{self.major_version_new}.{self.minor_version_new}" in self.versions:
            # create new version
            with self.connect() as swift:
                post = swift.post(container=f"{self.benchmark}.{self.major_version_new}.{self.minor_version_new}", options={'read_acl': ".r:*,.rlistings"})
                if not post['success']:
                    raise SwiftError(f"Version creation of {self.benchmark}.{self.major_version_new}.{self.minor_version_new} failed: {post['error']}")
                # except SwiftError as e:
                #     logger.error(e.value)
        
            # update
            self._get_versions()
        else:
            raise ValueError("Version already exists")

        if not f"{self.major_version_new}.{self.minor_version_new}" in self.versions:
            raise ValueError("Version creation failed")


    def _get_objects(self, readonly=False):
        if self.major_version is None or self.minor_version is None:
            raise ValueError("No version provided")
        if not 'preauthtoken' in self.auth_options.keys() or readonly:
            containername = f"{self.benchmark}.{self.major_version}.{self.minor_version}"
            url = f"{self.auth_options['preauthurl']}/{containername}"
            params={'format': 'xml'}
            response = requests.get(url, params=params)
            if response.ok:
                response_text = response.text
            else:
                response.raise_for_status()
            soup = BeautifulSoup(response_text, 'xml')
            names = soup.find_all('object')

            files = dict()
            for obj in names:
                mtime, accesstime = get_meta_mtime(self.auth_options['preauthurl'], containername, obj.find('name').text)
                files[obj.find('name').text] = {
                    "hash": obj.find('hash').text,
                    "size": obj.find('bytes').text,
                    "last_modified": obj.find('last_modified').text,
                    "symlink_path": obj.find('symlink_path').text if obj.find('symlink_path') is not None else "",
                    "x-object-meta-mtime": mtime,
                    "accesstime": accesstime
                }
            self.files = files
        
        else:
            with self.connect() as swift:
                files = dict()
                # get all objects
                for page in swift.list(container=f"{self.benchmark}.{self.major_version}.{self.minor_version}", options={'long' :  True}):
                    for element in page['listing']:
                        files[element['name']] = {'size': element['bytes'], 'last_modified': element['last_modified'], 'hash': element['hash']}
                        if 'symlink_path' in element.keys():
                            files[element['name']]['symlink_path'] = element['symlink_path']
                        else:
                            files[element['name']]['symlink_path'] = ""

                # stats for each object
                for element in swift.stat(container=f"{self.benchmark}.{self.major_version}.{self.minor_version}", objects=files.keys()):
                    # files[element['object']]['x-object-meta-last-modified'] = element['headers']['x-object-meta-last-modified']
                    files[element['object']]['x-object-meta-mtime'] = element['headers']['x-object-meta-mtime']

                self.files = files
        
    def find_objects_to_copy(self, reference_time = None, tagging_type = "all"):
            """
            Finds objects to copy based on the reference time and tagging type.

            Args:
                reference_time (datetime.datetime, optional): The reference time to compare with the object's time. Defaults to None.
                tagging_type (str, optional): The tagging type to consider. Defaults to "all".

            Raises:
                ValueError: If the tagging type is invalid or the reference time is not a datetime object.

            Returns:
                None
            """
            if tagging_type not in ["all"]:
                raise ValueError("Invalid tagging type")
            if tagging_type == "all":
                if reference_time is None:
                    reference_time = datetime.datetime.now()
                elif not type(reference_time) is datetime.datetime:
                    raise ValueError("Invalid reference time, must be datetime object")
                if len(self.files) == 0:
                    self._get_objects()

                if len(self.files) == 0:
                    return
                else:
                    for filename in self.files.keys():
                        if 'x-object-meta-mtime' in self.files[filename].keys():
                            file_time = datetime.datetime.fromtimestamp(float(self.files[filename]['x-object-meta-mtime']))
                        elif 'x-object-meta-last-modified' in self.files[filename].keys():
                            file_time = dateutil.parser.parse(self.files[filename]['x-object-meta-last-modified'])
                        else:
                            file_time = dateutil.parser.parse(self.files[filename]['last_modified'])
                        if file_time < reference_time:
                            self.files[filename]['copy'] = True
                        else:
                            self.files[filename]['copy'] = False
    
    def copy_objects(self, type = "copy"):
            """
            Copy objects from the current version to the new version.

            Args:
                type (str): The type of copying to perform. Valid values are "copy" and "symlink".
            """
            if type not in ["copy", "symlink"]:
                raise ValueError("Invalid type")
            if self.major_version is None or self.minor_version is None or self.major_version_new is None or self.minor_version_new is None:
                raise ValueError("No version provided")

            if type == "copy":
                with self.connect() as swift:
                    filenames = [filename for filename in self.files.keys() if self.files[filename]['copy']]
                    for out in swift.copy(f"{self.benchmark}.{self.major_version}.{self.minor_version}", filenames, {"destination": f"/{self.benchmark}.{self.major_version_new}.{self.minor_version_new}"}):
                        if out['success'] and out['action'] == 'copy_object':
                            self.files[out['object']]['copied'] = True
                        elif not out['success'] and out['action'] == 'copy_object':
                            self.files[out['object']]['copied'] = False
            if type == "symlink":
                raise NotImplementedError("Symlink copying not implemented")
    
    def create_new_version(self, version_type="minor", tagging_type="all", copy_type="copy"):
        """
        Wrapper to create new version and copy objects. 

        Args:
            version_type (str, optional): The type of version to create. Can be "minor" or "major". Defaults to "minor".
            tagging_type (str, optional): The type of tagging to apply. Defaults to "all".
            copy_type (str, optional): The type of copying to perform. Defaults to "copy".
        """
        if version_type not in ["minor", "major"]:
            raise ValueError("Invalid version type: 'minor' or 'major'")
        self.set_current_version()
        if version_type == "minor":
            self.set_new_version(None, True)
        elif version_type == "major":
            self.set_new_version(True, None)

        self._create_new_version()
        self._get_objects()
        self.find_objects_to_copy(tagging_type=tagging_type)
        self.copy_objects(copy_type)

    def archive_version(self, major_version, minor_version):
        NotImplementedError
        # self._update_overview(cleanup=True)

    def delete_version(self, major_version, minor_version):
        NotImplementedError
        # self._update_overview(cleanup=True)




