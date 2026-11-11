# -*- coding: utf-8 -*-
"""生成 A1 OSF 登记操作指引 Markdown（英文块与 a1_osf_fill_pack 同源）"""
import io
import os
import sys

HERE = os.path.dirname(__file__)
sys.path.insert(0, HERE)
from a1_osf_fill_pack import FIELDS  # noqa: E402

ROOT = os.path.join(HERE, '..', '..')
OUT = os.path.join(ROOT, '_markdown源', 'A1_OSF登记操作指引_20261111.md')

HEAD = """# A1 论文 OSF 回顾性登记操作指引

本件是操作手册，不是登记正文。**登记正文在 `A1_启动/results/a1_osf_paste_blocks.txt`**，按字段编号逐块复制粘贴；本件负责告诉你每一步点哪里、选什么、以及哪些动作不可逆。

全部内容基于 2026-11-11 实测：OSF 公开 schema API 取得 Secondary Data Preregistration 模板现行 27 个字段（schema 64775783798e08000a70407e，version 3，5 页），OSF 帮助中心核实提交、禁运、DOI 与撤回机制。

## 一、开做之前必须知道的三条规则

这三条决定操作方式，先看清再动手，否则会白做一次不可逆提交。

| 规则 | 具体后果 | 你的选择 |
| --- | --- | --- |
| 只有"立即公开"的注册才分配 DOI | 选禁运（embargo）的注册在禁运期内没有 DOI。而正文 Methods 与 Data availability 两处 [VERIFY] 需要填 OSF DOI，没有 DOI 就无法回填，登记也就失去意义 | 提交时选 Make registration public immediately，不勾 embargo |
| 提交后正文一字不可改 | OSF 只允许发 Update（须填理由），且 **Update 不能增删附件**。附件只能在提交时一次性带上 | 提交前逐字核对；附件在提交前全部上传 |
| 不要"从现有项目创建" | 该路径会把项目下所有文件（上限 5 GB）一并归档进注册，永久不可改 | 第一问选 No（从零创建），只附 4 个小文件 |

还有两条时间机制：提交后所有 admin 会收到确认邮件，**48 小时内无人操作则自动批准**；注册日期 = 提交日期，与批准日期无关。因此你单独一人操作时，提交即等于完成，不会卡在审批。

## 二、模板选择：为什么选 Secondary Data Preregistration

OSF 现行 45 个可用模板中，与本研究匹配的只有两个：

| 模板 | 是否选 | 理由 |
| --- | --- | --- |
| Secondary Data Preregistration | **选它** | 27 个字段分 5 页，其中第 4 页 "Knowledge of Data" 直接问 Prior Publication 与 Prior knowledge。这正是回顾性登记必须交代的两件事，模板给了正式位置 |
| OSF Preregistration（通用） | 不选 | 通篇按"数据尚未收集"设计，没有专门位置披露"我已经看过结果"，写在附加说明里会显得在遮掩 |
| Registered Report Protocol | 不可选 | 需先取得期刊 in-principle acceptance，本文尚在选刊阶段 |
| Prereg Challenge | 不可选 | 该模板已 inactive |

选对模板这件事本身就有评审价值：用二次数据模板做回顾性登记，等于在表单结构上就承认了数据是既有的。

## 三、分步操作

### Step 0　账号与身份（约 10 分钟，一次做完）

1. 用机构邮箱在 osf.io 注册，不要用 QQ/163 邮箱——OSF 的部分通知邮件会被国内邮箱判为垃圾邮件，而审批邮件一旦漏看会误事。
2. 进 Settings → Profile，做两件事：绑定 **ORCID**；填 **Affiliation**（广州医科大学附属中医医院）。这两项决定注册页的署名与机构检索归属，事后补要发 Update。
3. 记下你打算署名的英文形式。登记页与论文署名须一致，建议与投稿信一致：Yi Lu（鲁义主任医师）。

### Step 1　从零创建草稿（约 3 分钟）

1. 打开 osf.io/registries → 右上 **Add a Registration**。
2. 系统问 "Do you have content in an existing OSF project?" → 选 **No**。
3. 模板列表选 **Secondary Data Preregistration** → **Create Draft**。
4. OSF 会把草稿链接发到你的邮箱。**把邮件留着**，草稿也可在 osf.io/my-registrations?tab=drafts 找到。

草稿阶段可以反复改、可删除，不产生任何公开记录。真正不可逆的只有最后点 Register 那一下。

### Step 2　填元数据（约 15 分钟）

| 元数据项 | 建议填写 |
| --- | --- |
| Title | Power and phenotype-definition heterogeneity in drug-target Mendelian randomization: retrospective registration of an analysis plan for circulating CX3CL1 and diabetic neuropathy across five FinnGen releases |
| Description | 用下面这段，150 词以内，明确写 retrospective |
| Contributors | 先只加你自己并设为 Administrator。加合作者须其本人接受邀请，未接受会让审批卡住；发表前再加更稳 |
| License | 保留默认（CC BY 4.0）；若希望最大程度允许复用，选 CC0 1.0 |
| Subjects | Medical and Health Sciences → Epidemiology；另加 Genetics |
| Tags | Mendelian randomization; drug-target MR; diabetic neuropathy; CX3CL1; FinnGen; UK Biobank Pharma Proteomics Project; statistical power; negative results; reproducibility; STROBE-MR |

Description 建议文本：

ENBLOCK:
Retrospective registration of an analysis plan for a methodological study using
only publicly available summary-level data. This is not a prospective
preregistration. Exploratory analyses completed before registration are disclosed
in the registration and in any resulting manuscript. The study quantifies how much
conclusion strength, statistical power and verbal labels change across outcome
definitions, data releases and instrument sets in drug-target Mendelian
randomization, using circulating CX3CL1 and diabetic neuropathy across five
FinnGen releases as a worked example.
ENDENBLOCK

### Step 3　Study Information 页（字段 1—2）

两个长文本框。字段 2（Hypotheses）开头那段 "IMPORTANT FRAMING" 是全文最要紧的一段——它把"这些是结果已知后的预判"写死，是这次登记能否站住的关键。**一个字都不要删**。

### Step 4　Data Description 页（字段 3—11）

九项，其中三项是附件上传点（字段 9、11 与文件型字段）。三个要点：

- 字段 3 的样本量 **已核实为 n=33,436，全为欧洲裔**（Open Targets study 记录 `UKB_PPP_EUR_CX3CL1_P78423_OID20976_v1`，`nSamples=33,436`、`ldPopulation=nfe`、`relativeSampleSize=1.0`，2026-11-11 经 GraphQL API 核实，存档 `data/ot_cx3cl1_study_meta.json`）。该值与全部实算的 `N_EXP` 一致，直接粘贴即可，不再需要现场核对。
- 字段 6 的 URL 要特别注意 **R12 多一段 `/release/` 路径**，这是五个 release 里唯一的例外，抄错会让复现者直接 404。
- 字段 7 的取数日填 2026-11-10，与工程内所有取数文件的时间戳一致。

### Step 5　Variables 页（字段 12—18）

七项。字段 12（Manipulated variables）填 None 之后必须紧接着解释"暴露由遗传变异代理"，否则审稿人会以为填错模板。字段 15 必须写明三处缺失且不插补，这是可复核性的核心。字段 17 保留 LDLR 极端值做反例这一段，能在方法学审稿里加分。

### Step 6　Knowledge of Data 页（字段 19—20）

**这是全表最关键的一页，也是选这个模板的全部理由。**

字段 20（Prior knowledge）逐项列出登记前已完成的四类分析及其实测结果（首轮 MR P=0.764、clumping 后仅剩一个 cis、PP.H4=0.0584、六个 trans 位点）。不要缩写、不要含糊成"已做初步分析"——含糊比不写更糟，因为读者会假设最坏情况。

同页末尾主动写出三条"后果"（本登记不构成事前数据记录、方法学定位是事后决定的、登记后新增分析标为探索性），这三句是主动交代，不是被迫承认，写法上把它们放在预测性段落之前。

### Step 7　Analyses 页（字段 21—27）

七项，篇幅最大。三个必须写准的地方：

- 字段 21 主分析只写单工具变量 Wald ratio，其余全部列为敏感性分析；并写明 IVW 过原点、Egger 含截距。这两句是我们此前踩过坑的地方，写进去既准确也显示方法学自觉。
- 字段 23 开头必须写 "No a priori power or sample-size calculation was performed"。STROBE-MR 4b 就查这一句，把它写成"事后量化效能"而不是补救。
- 字段 24 的 Bonferroni 阈值写 0.05/133 = 3.76e-4，与四分类规则（I/II/III/IV）一并给出，并写明实测分布 0/0/15/4。这是本文的方法学产出，登记里给全，论文里就不必再解释规则来源。

### Step 8　提交设置（最容易出错的一步）

1. 点 **Review**，从头到尾读一遍。系统会把所有答案汇总在一页，这是最后一次改的机会。
2. 隐私选项选 **Make registration public immediately**。不要选 embargo。
3. 点 **Register** / **Submit**。
4. 立刻复制注册页 URL（形如 osf.io/ab3xy）存到本地。批准期 48 小时，期间页面显示 Pending Approval，仅贡献者可见。

### Step 9　取 DOI 并回填（约 15 分钟）

1. 公开后打开注册页，看 **Metadata** 区是否已显示完整 DOI（形如 10.17605/OSF.IO/AB3XY）。若显示的是 **Create DOI** 按钮，点它生成。部分注册需手动生成，别以为没有就是没分配。本注册实际 DOI 为 **10.17605/OSF.IO/APWE2**，位置在 overview 页最底部 Metadata 区的 License 下一行（Review 页不显示）。
2. 把 DOI 与注册日期回填到下列位置，逐项打勾。

## 四、提交后回填清单

| 序号 | 回填位置 | 填什么 | 状态 |
| --- | --- | --- | --- |
| 1 | 英文稿 Methods 段（M6） | registration date 与 DOI，替换原 [VERIFY] 占位 | 已完成：DOI 10.17605/OSF.IO/APWE2 已回填 |
| 2 | 英文稿 Declarations 段（DC） | OSF DOI，替换原 [VERIFY] 占位 | 已完成：DOI 10.17605/OSF.IO/APWE2 已回填 |
| 3 | `A1_启动/results/a1_verify_register.tsv` V02 行 | status 改 REGISTERED，note 填日期与 DOI | 已完成：DOI 10.17605/OSF.IO/APWE2 已回填 |
| 4 | `A1_启动/results/a1_version_snapshot.json` | 增加 osf_doi 与 osf_registered 字段 | 已完成：DOI 10.17605/OSF.IO/APWE2 已回填 |
| 5 | 投稿信 cover letter | 若提及预注册，补 DOI 与"回顾性"字样 | 已完成：DOI 10.17605/OSF.IO/APWE2 已回填 |
| 6 | 选刊确定后 | 三方互链：OSF ↔ GitHub 仓库 ↔ Zenodo DOI | 待做（需 GitHub/Zenodo 账号动作） |

## 五、27 个字段逐条

每一条给中文提示与可直接粘贴的英文正文。编号与 OSF 表单从上到下的顺序一致。完整纯文本版另存 `a1_osf_paste_blocks.txt`，直接从 txt 复制比从 Word 复制更安全（不会带入引号与样式）。

"""

TAIL = """
## 六、附件策略

注册提交后**不能增删附件**，所以现在就要定。只附四份体积小、指向明确的文件，合计不到 200 KB：

| 附件 | 作用 | 是否必附 |
| --- | --- | --- |
| a1_endpoint_dict.tsv | 19 个组合的端点字典，字段 3、6、10 指向它 | 必附 |
| a1_fetch_log.json | 133 次取数的逐请求状态与缺失变异记录，证明取数可复核 | 必附 |
| a1_tabix_fetch.py | 定点取数实现，字段 9 的文档 | 必附 |
| a1_multi_fetch.py | 批量取数与 harmonise 实现，字段 9 的文档 | 建议附 |

不要附原始 GWAS 汇总统计、不要附整包代码、不要附论文全文。原始数据在公开桶里，代码将来走 GitHub + Zenodo，OSF 只承担时间戳与登记功能。附得越多，日后越难把"登记"和"发表"分开。

## 七、提交前核对 12 项

打印这一页，逐项打勾。任何一项不确定，就不要点 Register。

| 核对项 | 打勾 |
| --- | --- |
| 模板确为 Secondary Data Preregistration（页数与字段数对得上：5 页 27 字段） | □ |
| 第一问选了 No，没有把任何现有项目的文件带进注册 | □ |
| 标题无拼写错误，且含 retrospective registration 字样 | □ |
| 27 个字段全部有内容，无红框或感叹号 | □ |
| 字段 2 首段 IMPORTANT FRAMING 完整保留 | □ |
| 字段 20 Prior knowledge 四项分析齐全，数字与内部报告一致 | □ |
| 字段 6 抽查三个 URL 能在浏览器打开（含 R12 的 /release/ 路径） | □ |
| 字段 7 取数日填 2026-11-10 | □ |
| 字段 3 的样本量已现场核实，或已改为文字表述 | □ |
| 隐私选项为立即公开，未勾 embargo | □ |
| 四份附件已上传完毕 | □ |
| 草稿链接与注册页 URL 已存本地 | □ |

## 八、红线

以下四件事在任何情况下都不做，做了这次登记反而变成负资产：

- 删掉 IMPORTANT FRAMING 或 Prior knowledge 两段，把回顾性写成前瞻性。一旦被核查时间戳与内部报告日期，即为 HARKing 实证，比不做登记更严重。
- 为了让登记"看起来干净"而调整任何 P 值、OR 或置信区间。登记里的数字必须与结果文件逐位一致。
- 提交后直接新建第二个注册来"修正"。任何偏离走 Update 流程，附理由，保留修订痕迹。
- 在登记里写入尚未核实的样本量、release 官方公开日期或期刊政策。核不到就写文字表述或留 [VERIFY]，让它留白。

## 九、本轮产出

| 文件 | 内容 |
| --- | --- |
| results/a1_osf_paste_blocks.txt | 27 个字段的纯文本粘贴包，本件的操作对象 |
| results/a1_osf_field_mapping.tsv | 字段名、输入类型、中文提示、英文值、字数，可复核 |
| results/a1_osf_registration_pack.txt | 此前备好的总纲式登记文本（保留，用于投稿信与论文 Methods 引用） |
| data/osf_schema_secondary.json | OSF 模板 schema 原始快照，证明字段名与顺序取自现行版本 |

口径声明：登记内容只陈述已完成的分析与已核实的数值，不预测结果、不承诺结论方向，不改变"未获得支持性证据、现有公开数据不足以检验中等效应"这一表述，不推翻标书假说。登记是合规与透明的底线动作，不构成对研究结论的背书。
"""


def main():
    parts = [HEAD]
    page_now = None
    for i, d in enumerate(FIELDS, 1):
        if d['page'] != page_now:
            page_now = d['page']
            parts.append('### %s 页\n\n' % page_now)
        parts.append('#### %02d　%s\n\n' % (i, d['field']))
        parts.append('%s\n\n' % d['cn_hint'])
        if d['input_type'] == 'file upload':
            parts.append('（本字段不上传正文，按第六节附件策略处理。）\n\n')
            continue
        parts.append('ENBLOCK:\n%s\nENDENBLOCK\n\n' % d['en'].strip())
    parts.append(TAIL)
    md = ''.join(parts)
    with io.open(OUT, 'w', encoding='utf-8') as fh:
        fh.write(md)
    print('written', OUT)
    print('chars', len(md), 'enzblocks', md.count('ENBLOCK:'))


if __name__ == '__main__':
    main()
