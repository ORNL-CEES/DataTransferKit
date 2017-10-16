# Building DTK
- To build DTK on ALCF ressources see [README_ALCF.md](README_ALCF.md)
- To build DTK on Summitdev see [README_OLCF_SUMMITDEV.md](README_OLCF_SUMMITDEV.md)
- To build DTK on CADES ressources see [README_CADES.md](README_CADES.md)

# Changelog generation
We are using [github_changelog_generator](https://github.com/skywinder/github-changelog-generator) which is a Ruby-based fully automatic changelog generator based on **tags**, **issues** and merged **pull requests**.

## Installation
Install `rubygems` and run
```bash
$ gem install github_changelog_generator
```
If your Ruby is too old, you can use a docker image
```bash
$ docker run -it --rm \
        -v $(pwd):/app \
        prooph/github-changelog-generator \
        ORNL-CEES/DataTransferKit \
        --no-pull-requests \
        --include-labels 'bug,enhancement,New Feature' \
        --enhancement-labels 'enhancement,New Feature'
```
