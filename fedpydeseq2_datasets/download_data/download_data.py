import argparse
import shutil
import subprocess
import tempfile
from pathlib import Path

import yaml  # type: ignore


def download_data(
    config_path: str | Path,
    download_data_directory: str | Path,
    raw_data_output_path: str | Path,
    snakemake_env_name: str,
    conda_activate_path: str | Path | None = None,
):
    """
    Download the data.

    Parameters
    ----------
    config_path : Union[str, Path]
        The path to the configuration file.
    download_data_directory : Union[str, Path]
        The path to the download data directory.
    raw_data_output_path : Union[str, Path]
        The path to the raw data output directory.
    snakemake_env_name : str
        The name of the Snakemake environment.
    conda_activate_path : Optional[Union[str, Path]], optional
        The path to the conda activate script, by default None.
    """
    temp_dir = tempfile.mkdtemp()
    try:
        # Copy the specified directory to the temporary directory
        tmp_download_data_path = Path(temp_dir, "download_data")
        shutil.copytree(download_data_directory, tmp_download_data_path)

        # Open the configuration file
        with open(config_path) as f:
            config = yaml.safe_load(f)

        # Add a field "output_path" to the configuration with the raw data output path
        config["output_path"] = str(raw_data_output_path)
        config_file = Path(tmp_download_data_path, "config", "config.yaml")
        config_file.parent.mkdir(parents=True, exist_ok=True)
        with open(config_file, "w") as f:
            yaml.dump(config, f)

        print("Config file", config_file)

        # Create a conda env
        # make an envs direcrory in the temporary directory
        envs_dir = Path(temp_dir, "envs")
        envs_dir.mkdir(parents=True, exist_ok=True)

        env_prefix = Path(envs_dir, snakemake_env_name)

        create_conda_env(
            env_file=Path(tmp_download_data_path, "snakemake_env.yaml"),
            env_prefix=env_prefix,
        )

        access_conda_command = (
            f"""
        cd {tmp_download_data_path}
        echo `{tmp_download_data_path}`
        conda init bash
        if [ -f ~/.bashrc ]; then
            . ~/.bashrc
        fi
        if [ -f ~/.bash_profile ]; then
            . ~/.bash_profile
        fi
        """
            if conda_activate_path is None
            else f"""
        . {conda_activate_path}
        """
        )
        command = f"""
        {access_conda_command}
        conda activate {env_prefix}
        cd {tmp_download_data_path}
        snakemake --cores all --resources download_concurrent=3 \
              --use-conda
        """

        subprocess.run(command, shell=True, check=True, executable="/bin/bash")

    except Exception as e:  # noqa BLE001
        print(f"An error occurred: {e}")
    finally:
        # Clean up the temporary directory
        shutil.rmtree(temp_dir)


def create_conda_env(env_file: str | Path, env_prefix: str | Path):
    """
    Create a Conda environment.

    Parameters
    ----------
    env_file : str or Path
        The path to the environment file.
    env_prefix : str or Path
        The prefix (location) where the Conda environment will be created.
    """
    try:
        # Create the Conda environment
        create_env_cmd = (
            f"yes | conda env create --prefix {env_prefix} --file {env_file}"
        )
        subprocess.run(create_env_cmd, shell=True, check=True, executable="/bin/bash")
        print(f"Conda environment created successfully at '{env_prefix}'.")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while creating the Conda environment: {e}")


def main():
    """Run main function to download the data."""
    parser = argparse.ArgumentParser("""Download the data.""")
    parser.add_argument(
        "--only_luad", action="store_true", help="Only download the LUAD dataset"
    )
    parser.add_argument(
        "--raw_data_output_path",
        type=str,
        required=False,
        help="Path to the raw data output directory",
    )
    parser.add_argument(
        "--conda_activate_path",
        type=str,
        help="Path to the conda activate script",
        required=False,
    )
    args = parser.parse_args()

    if args.only_luad:
        config_path = Path(__file__).parent / "config" / "config_luad.yaml"
    else:
        config_path = Path(__file__).parent / "config" / "config.yaml"

    if args.conda_activate_path is not None:
        conda_activate_path = Path(args.conda_activate_path)
    else:
        conda_activate_path = None

    if args.raw_data_output_path is None:
        raw_data_output_path = Path(__file__).parent.parent.parent / "data/raw"
    else:
        raw_data_output_path = Path(args.raw_data_output_path)
    download_data_directory = Path(__file__).parent

    download_data(
        config_path=config_path,
        download_data_directory=download_data_directory,
        raw_data_output_path=raw_data_output_path,
        snakemake_env_name="snakemake_env",
        conda_activate_path=conda_activate_path,
    )


if __name__ == "__main__":
    main()
