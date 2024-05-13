<h1 align="center">Leash Bio Kaggle</h1>

<h4 align="center">Predict small molecule-protein interactions using the Big Encoded Library for Chemical Assessment (BELKA)</h4>

## Deploying

We use [bump-my-version](https://github.com/callowayproject/bump-my-version) to release a new version.
This will create a git tag that is used by [poetry-dynamic-version](https://github.com/mtkennerly/poetry-dynamic-versioning) to generate version strings and update `CHANGELOG.md`.

For example, to bump the `minor` version you would run the following command.

```bash
poetry run bump-my-version bump minor
```

After releasing a new version, you need to push and include all tags.

```bash
git push --follow-tags
```

## License

This project is released under the MIT License as specified in `LICENSE.md`.
