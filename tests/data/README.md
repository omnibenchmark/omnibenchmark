# tests data

## bundles

To avoid cloning over the network, and to achieve more self-contained tests, we have a few git bundles embedded in this repo.

There's a script to create/update bundles from existing network:

```
python bundle_repos.py
```

When creating and committing bundles, please take into account that there's a
risk of adding too much cruft to this repo history. So practice
self-moderation! Test repos should be as slim as possible.

Also, avoid commiting the repos folder.
