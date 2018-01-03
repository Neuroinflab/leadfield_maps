import parameters as params
import requests
import tarfile
import os

#Low res
points_url = 'https://www.dropbox.com/s/bmjc08m9dgi9vha/points.tar.gz?dl=1'
sigmas_url = 'https://www.dropbox.com/s/8vnrh66m6goh409/sigmas.tar.gz?dl=1'
mesh_url = 'https://www.dropbox.com/s/413x4dhxy7g9fle/mesh.tar.gz?dl=1'


def extract_url(url, name, path):
    r = requests.get(url, allow_redirects=True)
    fname = name + '.tar.gz'
    open(fname, 'wb').write(r.content)
    tar = tarfile.open(fname, "r:gz")
    tar.extractall(path)
    tar.close()
    os.remove(fname)
    print('Done fetching for: ', name)
    return


extract_url(points_url, 'points', params.points_path)
extract_url(sigmas_url, 'sigmas', params.sigmas_path)
extract_url(mesh_url, 'mesh', params.mesh_path)
