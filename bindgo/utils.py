import numpy as np
import numpy.typing as npt


def split_indices(
    n_rows: int, train_split: float = 0.8, size: int | None = None
) -> tuple[npt.NDArray[np.uint64], npt.NDArray[np.uint64]]:
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
    indices = np.arange(n_rows, dtype=np.uint64)
    np.random.shuffle(indices)

    # Trim the number of indices to size.
    indices = indices[None:size]

    # Split indices into training and validation sets
    train_size = int(indices.shape[0] * train_split)
    train_indices = indices[:train_size]
    val_indices = indices[train_size:]

    return train_indices, val_indices


def data_generator(scanner, indices, batch_size=1000):
    indices_set = set(indices)
    for i, batch in enumerate(scanner.to_batches()):
        start_idx = i * batch.num_rows
        end_idx = start_idx + batch.num_rows
        batch_indices = range(start_idx, end_idx)
        selected_indices = [idx for idx in batch_indices if idx in indices_set]
        if selected_indices:
            selected_rows = batch.take(selected_indices)
            yield selected_rows
