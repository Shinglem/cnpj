import com.github.jengelman.gradle.plugins.shadow.tasks.ShadowJar
import org.jetbrains.kotlin.allopen.gradle.AllOpenExtension
import org.jetbrains.kotlin.gradle.tasks.KotlinCompile
import org.jetbrains.kotlin.noarg.gradle.NoArgExtension

plugins {
    java
    kotlin("jvm")
    id("org.jetbrains.kotlin.plugin.noarg")
    id("org.jetbrains.kotlin.plugin.allopen")
    id("com.github.johnrengelman.shadow")
}

group = "io.github.shinglem"

repositories {
    mavenCentral()
}
val vertxVersion: String by project
val springBootVersion: String by project
dependencies {

    implementation(group = "io.vertx", name = "vertx-web", version = vertxVersion)
    implementation(group = "io.vertx", name = "vertx-core", version = vertxVersion)
    implementation(group = "io.vertx", name = "vertx-lang-kotlin", version = vertxVersion)
    implementation(group = "io.vertx", name = "vertx-lang-kotlin-coroutines", version = vertxVersion)

    implementation(group= "ch.qos.logback", name= "logback-core", version= "1.2.3")
    implementation(group= "ch.qos.logback", name= "logback-classic", version= "1.2.3")
    implementation(group= "org.slf4j", name= "slf4j-api", version= "1.7.28")

    implementation(group = "org.jetbrains.kotlinx", name = "kotlinx-coroutines-jdk8", version = "1.2.2")
    implementation(kotlin("stdlib-jdk8"))
    testCompile("junit", "junit", "4.12")
}

configure<JavaPluginConvention> {
    sourceCompatibility = JavaVersion.VERSION_1_8
}
tasks.withType<KotlinCompile> {
    kotlinOptions.jvmTarget = "1.8"
}


configure<AllOpenExtension> {
    this.annotations(
        "org.springframework.transaction.annotation.Transactional",
        "org.springframework.stereotype.Component",
        "org.springframework.stereotype.Controller",
        "org.springframework.stereotype.Indexed",
        "org.springframework.stereotype.Repository",
        "org.springframework.stereotype.Service"
    )
}

configure<NoArgExtension> {
    this.annotations(
        "com.sanss.ibp.bac.common.util.NoArg"
    )
}

tasks.withType<ShadowJar>{
    archiveBaseName.set("shadow")
    archiveClassifier.set("")
    archiveVersion.set("")
    manifest {
        this.attributes ("Main-Class" to "io.github.shinglem.cnpj.AppKt")
    }
}