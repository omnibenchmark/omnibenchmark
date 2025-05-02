# Documentation

This repo serves: [https://omnibenchmark.org](https://omnibenchmark.org)

## GitLab CI/CD

This project's static Pages are built by [GitLab CI/CD](https://about.gitlab.com/stages-devops-lifecycle/continuous-integration/),
following the steps defined in [`.gitlab-ci.yml`](.gitlab-ci.yml).

The [`.gitlab-ci.yml`](.gitlab-ci.yml) has a lengthy documentation discussing the automated deployment (staging) and the manual deployment to production. In short:
- edits in non master/main branchs will trigger a deployment to `review.omnibenchmark.org`, which can be aborted
- edits in master/main branches will trigger a deployment to `staging.omnibenchmark.org`, which can be manually confirmed to be moved to `production.omnibenchmark.org` and its synonym `omnibenchmark.org`


## Building locally

To work locally with this project (which might not be necessary, as the gitlab CI/CD is already designed so it won't deploy to production but to review.omnibenchmark.org or to staging.omnibenchmark.org instead), you'll have to follow the steps below:

1) Fork, clone or download this project.

2) Install the Python modules in `requirements.txt`.

To install them inside a virtenv you can:

```
## check your path to deactivate conda envs or anything weird there
echo $PATH

# once the path is clean (i.e. with conda deactivate or what needed)
# create a folder to store virtenvs, if not existing, and go there
mkdir -p ~/virtenvs
cd $_

# create a python3 virtualenv with your system's python3
python3 -m venv mkdocs

# activate the virtualenv
source ~/virtenvs/mkdocs/bin/activate

# double check: which pip are we using now?
# should be the one within the `mkdocs` virtenv
which pip 

# go to your home to clone the omni documentation repo there
cd

# clone it (master/main branch)
git clone git@gitlab.renkulab.io:omnibenchmark/omni_site.git

# go to the folder, install using the requirements.txt file
cd omni_site
pip install -r requirements.txt
mkdocs --version

# deactivate your virtenv when finished
deactivate
```

3) Clone the repository and bring your modifications to files in `docs`. If you wish to rearrange or reorder pages, beware to also change `mkdocs.yml`. This file defines how pages are ordered and which files they are based on.

4) When you are done, **test** your modifications locally by running: 

    `mkdocs serve` 

5) commit and push your changes (but **not the HTMLs**). 


Some useful links: 

- [MkDocs](https://www.mkdocs.org/): basic documentation for MkDocs.

- [MkDocs material](https://squidfunk.github.io/mkdocs-material/): documentation for advanced customization, mainly `themes` and `features` from `mkdocs.yml`. 

- [MkDocs - RealPython](https://realpython.com/python-project-documentation-with-mkdocs/#step-6-host-your-documentation-on-github): explains how to use github for documentation and how to deploy it. 
