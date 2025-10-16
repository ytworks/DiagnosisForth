#!/usr/bin/env python3
"""
テスト用ログ送信スクリプト
既存のlogger.pyの関数を使用してSlackにテストログを送信します
"""
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from DiagnosisForth.utils.logger import log_sdf, log_pdb, log_pharmacophore, _send_log


def test_all_logging_functions():
    """全てのロギング関数をテストする"""
    print("=== ログ送信テスト開始 ===\n")
    
    # 1. 直接メッセージ送信テスト
    print("1. 直接メッセージ送信テスト...")
    _send_log("テストメッセージ: logger.pyの動作確認")
    print("   ✓ 送信完了")
    time.sleep(1)
    
    # 2. log_sdf テスト
    print("\n2. log_sdf テスト...")
    log_sdf(
        file_path="/test/path/sample_molecule.sdf",
        smiles="CC(C)C[C@H](NC(=O)[C@H](Cc1ccccc1)NC(=O)OC(C)(C)C)C(=O)O",
        function_name="test_log_sdf"
    )
    print("   ✓ SDFログ送信完了")
    time.sleep(1)
    
    # 3. log_pdb テスト
    print("\n3. log_pdb テスト...")
    log_pdb(
        file_path="/test/path/protein_structure.pdb",
        function_name="test_log_pdb"
    )
    print("   ✓ PDBログ送信完了")
    time.sleep(1)
    
    # 4. log_pharmacophore テスト
    print("\n4. log_pharmacophore テスト...")
    log_pharmacophore(
        smiles="CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        function_name="test_log_pharmacophore"
    )
    print("   ✓ Pharmacophoreログ送信完了")
    
    print("\n=== 全テスト完了 ===")
    print("Slackのチャンネルを確認してください")


def test_custom_message():
    """カスタムメッセージを送信する"""
    import argparse
    parser = argparse.ArgumentParser(description="カスタムメッセージをSlackに送信")
    parser.add_argument("--message", "-m", type=str, help="送信するメッセージ")
    args = parser.parse_args()
    
    if args.message:
        print(f"カスタムメッセージ送信中: {args.message}")
        _send_log(args.message)
        print("✓ 送信完了")
    else:
        test_all_logging_functions()


if __name__ == "__main__":
    test_custom_message()