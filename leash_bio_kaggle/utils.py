import numpy as np
import numpy.typing as npt


def split_indices(
    n_rows, train_split: float = 0.8, size: str | None = None
) -> (npt.NDArray[np.uint64], npt.NDArray[np.uint64]):
    """
    Splits the indices of rows into training and validation sets.

    Args:
        n_rows: The total number of rows to generate indices for.
        train_split: The proportion of indices to allocate to the training set.
        size: If provided, trims the number of indices to this size.

    Returns:
        A tuple containing two arrays:

            - train_indices: Indices for the training set.
            - val_indices: Indices for the validation set.
    """

    # Generate indices and shuffle them in place.
    indices = np.arange(n_rows)
    np.random.shuffle(indices)

    # Trim the number of indices to size.
    indices = indices[None:size]

    # Split indices into training and validation sets
    train_size = int(indices.shape[0] * train_split)
    train_indices = indices[:train_size]
    val_indices = indices[train_size:]

    return train_indices, val_indices
