# 使用 Ubuntu 作为基础镜像
FROM ubuntu:24.04

# 作者和版本信息
LABEL Author="sherlock" \
      Version="v1.0"

# 设置环境变量
ENV PATH="/usr/local/bin:/opt/python-venv/bin:/opt/FEELnc/bin:/opt/FEELnc/scripts:$PATH" \
    R_LIBS_USER="/usr/local/lib/R/site-library" \
    VIRTUAL_ENV="/opt/python-venv" \
    PERL5LIB="$PERL5LIB:/opt/FEELnc/lib/" \
    DATA_DIR="/opt/pipeline/data"

# 复制数据目录和主脚本
COPY docker/data /opt/pipeline/data
COPY docker/main.sh /opt/pipeline

# 创建符号链接以便于执行
RUN ln -s /opt/pipeline/main.sh /usr/local/bin/lnc_pipeline

# 安装系统依赖项和必需工具，并清理缓存
RUN apt-get update && apt-get install -y \
    build-essential \
    zlib1g-dev \
    libncurses-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    wget \
    curl \
    ca-certificates \
    git \
    python3 \
    python3-venv \
    perl \
    bioperl \
    cpanminus \
    r-base \
    r-base-dev \
    r-cran-tidyverse \
    r-cran-stringr \
    r-cran-optparse \
    samtools \
    hisat2 \
    stringtie \
    gffread \
    cmake  # 添加 cmake 以编译 KmerInShort

RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# 安装 Python 2.7 并清理中间文件
RUN cd /usr/src && \
    wget https://www.python.org/ftp/python/2.7.18/Python-2.7.18.tgz && \
    tar xzf Python-2.7.18.tgz && \
    cd Python-2.7.18 && \
    ./configure --enable-optimizations && \
    make altinstall && \
    cd .. && \
    rm -rf Python-2.7.18 Python-2.7.18.tgz

# 安装 pip for Python 2.7
RUN curl https://bootstrap.pypa.io/pip/2.7/get-pip.py -o get-pip.py && \
    /usr/local/bin/python2.7 get-pip.py && \
    rm get-pip.py

# 设置 Python 3 虚拟环境并安装软件包
RUN python3 -m venv /opt/python-venv && \
    . /opt/python-venv/bin/activate && \
    pip install --upgrade pip && \
    pip install numpy pandas

# 在 Python 2.7 环境中安装 CPAT 和依赖
RUN /usr/local/bin/python2.7 -m pip install numpy pandas && \
    /usr/local/bin/python2.7 -m pip install cpat==1.2.4

# 安装 Perl 模块
RUN cpanm --force Bio::Perl && \
    cpanm --force Bio::DB::SeqFeature::Store && \
    cpanm Statistics::Basic && \
    cpanm List::MoreUtils && \
    cpanm Parallel::ForkManager

# 安装 KmerInShort
RUN git clone --recursive https://github.com/rizkg/KmerInShort /opt/KmerInShort && \
    cd /opt/KmerInShort && mkdir build && cd build && \
    cmake .. && make -j 8 && \
    chmod +x /opt/KmerInShort/build/KmerInShort && \
    ln -s /opt/KmerInShort/build/KmerInShort /usr/local/bin/KmerInShort

# 安装 fasta_ushuffle
RUN git clone https://github.com/agordon/fasta_ushuffle.git /opt/fasta_ushuffle && \
    cd /opt/fasta_ushuffle && make && \
    chmod +x /opt/fasta_ushuffle/fasta_ushuffle /opt/fasta_ushuffle/ushuffle && \
    ln -s /opt/fasta_ushuffle/fasta_ushuffle /usr/local/bin/fasta_ushuffle && \
    ln -s /opt/fasta_ushuffle/ushuffle /usr/local/bin/ushuffle

# 安装 FEELnc
RUN git clone https://github.com/tderrien/FEELnc.git /opt/FEELnc && \
    chmod +x /opt/FEELnc/bin/* /opt/FEELnc/scripts/*

# 安装 R 软件包并设置默认镜像
RUN Rscript -e 'options(repos = "http://cran.r-project.org"); install.packages(c("LncFinder", "seqinr", "optparse", "data.table", "ROCR", "randomForest"))'

# 安装 Diamond 工具
RUN wget https://github.com/bbuchfink/diamond/releases/download/v2.1.10/diamond-linux64.tar.gz -O diamond-linux64.tar.gz && \
    tar -xzf diamond-linux64.tar.gz && \
    mv diamond /usr/local/bin/ && \
    rm diamond-linux64.tar.gz

# 清理 APT 缓存和不必要的文件以减少镜像大小
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /usr/src/Python-2.7.18

# 定义默认命令
CMD ["lnc_pipeline"]
